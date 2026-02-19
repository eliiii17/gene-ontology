from abc import ABC, abstractmethod

import numpy as np


class SimilarityStrategy(ABC):
    """The interface for similarity algorithms between terms."""

    @abstractmethod
    def calculate_similarity(self, term_id_a, term_id_b, graph):
        """Calculates a similarity score and returns it(float)."""
        pass


class JaccardStrategy(SimilarityStrategy):
    """
    Computes similarity based on graph structure overlap.

    Jaccard(A, B) =
    |intersect(A.ancestors, B.ancestors)| / |union(A.ancestors, B.ancestors)|

    0.0 = No overlap, 1.0 = Identical terms.
    """

    def calculate_similarity(self, term_id_a, term_id_b, graph):
        ancestors_a = graph.get_ancestors(term_id_a)
        ancestors_b = graph.get_ancestors(term_id_b)

        intersects_num = len(ancestors_a.intersection(ancestors_b))
        unions_num = len(ancestors_a.union(ancestors_b))

        if unions_num == 0:
            return 0

        return intersects_num / unions_num


class WuPalmerStrategy(SimilarityStrategy):
    """
    Computes similarity based on the depth of the lowest(deepest) common ancestor(LCA).

    i.e. nodes sharing a deep ancestor make them very similar.

    WUP(A, B) = (2 * depth(LCA)) / (depth(A) + depth(B))
    """

    def calculate_similarity(self, term_id_a, term_id_b, graph):
        ancestors_a = graph.get_ancestors(term_id_a)
        ancestors_b = graph.get_ancestors(term_id_b)

        # Find common ancestors and the lowest/deepest one
        commons = ancestors_a.intersection(ancestors_b)
        if not commons:
            return 0.0

        # Calculate depths (distance from root)
        lca_depth = max(graph.get_depth(node) for node in commons)
        depth_a = graph.get_depth(term_id_a)
        depth_b = graph.get_depth(term_id_b)

        if (depth_a + depth_b) == 0:
            return 0.0

        return (2.0 * lca_depth) / (depth_a + depth_b)


class InformationContentCalculator:
    """
    Statistical model for information content.
    Calculates how informative(specific or rare) each GO term is.
    """

    def __init__(self, graph, repo):
        # this is a dict, which maps each GO term ID to its information content value
        self._ic_map = {}
        self._compute_ic(graph, repo)

    def _compute_ic(self, graph, repo):
        """Builds the information content (IC) map using numpy."""

        # count total number of genes
        N_GENES = len(repo.genes)

        # count genes for every term, and store it here
        # dict like GO id: number of related genes
        term_counts = {}

        for term_id in graph.all_term_ids():
            genes = repo.get_genes_for_term_recursive(term_id, graph)

            if len(genes) > 0:
                term_counts[term_id] = len(genes)

        # convert to numpy
        term_counts_array = np.array(list(term_counts.values()), dtype=float)

        # avoiding log(0) issues
        term_counts_array[term_counts_array == 0] = 1

        # probabilities: p(term) = count / total
        probs = term_counts_array / N_GENES

        # information content (ic) = -log(p)
        ic_values = -np.log(probs)

        # get IDs back
        keys = list(term_counts.keys())
        self._ic_map = dict(zip(keys, ic_values))

    def get_ic(self, term_id):
        return self._ic_map[term_id]


class ResnikStrategy(SimilarityStrategy):
    """
    Semantic similarity based on information content.
    Requires an InformationContentCalculator instance.
    """

    def __init__(self, ic_calculator):
        self._calculator = ic_calculator

    def calculate_similarity(self, term_id_a, term_id_b, graph):
        """
        In Resnik strategy, the similarity score is the maximum information content
        of the common ancestors
        """
        ancestors_a = graph.get_ancestors(term_id_a)
        ancestors_b = graph.get_ancestors(term_id_b)

        common_ancestors = ancestors_a.intersection(ancestors_b)

        if not common_ancestors:
            return 0

        max_ic = 0

        for ancestor in common_ancestors:
            ic = self._calculator.get_ic(ancestor)

            if ic > max_ic:
                max_ic = ic

        return max_ic


class GeneSimilarityCalculator:
    """
    Calculates similarity between two genes based on their GO annotations
    by pre-computing term-pair scores.
    """

    def __init__(self, term_strategy, graph):
        self._strategy = term_strategy
        self._graph = graph

    def calculate_similarity(self, gene_a, gene_b):
        """
        Calculates the similarity between two genes.
        Returns a tuple of (overall_score, top_5_term_matches).
        """
        terms_a = gene_a.get_valid_go_ids()
        terms_b = gene_b.get_valid_go_ids()

        if not terms_a or not terms_b:
            return (0.0, [])

        scores = self._precompute_scores(terms_a, terms_b)

        score_ab = self._avg_best_match(terms_a, terms_b, scores)
        score_ba = self._avg_best_match(terms_b, terms_a, scores)
        overall_score = (score_ab + score_ba) / 2.0

        details = self._get_top_details(scores)

        return (overall_score, details)

    def _get_pair_key(self, t1, t2):
        return tuple(sorted((t1, t2)))

    def _precompute_scores(self, terms_a, terms_b):
        scores = {}
        for ta in terms_a:
            for tb in terms_b:
                key = self._get_pair_key(ta, tb)
                if key not in scores:
                    scores[key] = self._strategy.calculate_similarity(
                        ta, tb, self._graph
                    )
        return scores

    def _avg_best_match(self, source_terms, target_terms, scores):
        if not source_terms:
            return 0.0

        best_scores = []
        for s in source_terms:
            # Find best similarity for term s against any term in target_terms
            s_scores = [scores.get(self._get_pair_key(s, t), 0.0) for t in target_terms]
            best_scores.append(max(s_scores) if s_scores else 0.0)

        return sum(best_scores) / len(source_terms)

    def _get_top_details(self, scores):
        # Create list of (term_a, term_b, score) for non-identical pairs
        valid_pairs = [
            (t1, t2, score) for (t1, t2), score in scores.items() if t1 != t2
        ]
        valid_pairs.sort(key=lambda x: x[2], reverse=True)

        details = []
        for t1, t2, score in valid_pairs[:5]:
            term_a = self._graph.get_term(t1)
            term_b = self._graph.get_term(t2)
            details.append(
                {
                    "term_a": t1,
                    "term_a_name": term_a.name if term_a else "N/A",
                    "term_b": t2,
                    "term_b_name": term_b.name if term_b else "N/A",
                    "score": score,
                }
            )
        return details

    def calculate_matrix(self, genes):
        """
        Calculates a similarity matrix for a list of genes.
        Returns a tuple (labels, matrix) where:
        labels is a list of gene symbols
        matrix is a numpy array (N, N). N is the number of genes
        """
        n = len(genes)
        matrix = np.zeros((n, n))
        labels = [gene.symbol for gene in genes]

        # Fill the matrix
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    score = 1.0  # self similarity set to 1
                else:
                    score, _ = self.calculate_similarity(genes[i], genes[j])

                matrix[i, j] = score
                matrix[j, i] = score

        return labels, matrix


class SimilarityService:
    """
    Handles similarity calculations.
    Encapsulates strategies and handles request processing.
    """

    def __init__(self, graph, repo, ic_calculator):
        self.graph = graph
        self.repo = repo

        # Initialize strategies
        self.strategies = {
            "jaccard": JaccardStrategy(),
            "wupalmer": WuPalmerStrategy(),
            "resnik": ResnikStrategy(ic_calculator),
        }

        # Initialize gene calculators for each strategy
        self.gene_calculators = {
            name: GeneSimilarityCalculator(strategy, graph)
            for name, strategy in self.strategies.items()
        }

    def process_request(self, args):
        """
        Processes a similarity request (from Flask args) and returns the view data.
        """
        mode = args.get("mode", "term")
        strategy_name = args.get("strategy", "jaccard")

        term_a_query = args.get("term_a_query")
        term_b_query = args.get("term_b_query")
        term_a_id = args.get("term_a_id")
        term_b_id = args.get("term_b_id")
        gene_a_query = args.get("gene_a_query")
        gene_b_query = args.get("gene_b_query")
        gene_list_query = args.get("gene_list_query")

        view_data = {
            "mode": mode,
            "strategy": strategy_name,
            "term_a_query": term_a_query,
            "term_b_query": term_b_query,
            "term_a_id": term_a_id,
            "term_b_id": term_b_id,
            "gene_a_query": gene_a_query,
            "gene_b_query": gene_b_query,
            "gene_list_query": gene_list_query,
            "result": None,
            "error": None,
            "details": None,
            "matrix_data": None,
        }

        # Validate Strategy
        if strategy_name not in self.strategies:
            view_data["error"] = f"Unknown strategy: {strategy_name}"
            return view_data

        if mode == "term" and term_a_id and term_b_id:
            term_a = self.graph.get_term(term_a_id)
            term_b = self.graph.get_term(term_b_id)

            if not term_a or not term_b:
                view_data["error"] = (
                    "One or both GO terms not found. Please select a valid term from the list."
                )
            else:
                strategy = self.strategies[strategy_name]
                view_data["result"] = strategy.calculate_similarity(
                    term_a.id, term_b.id, self.graph
                )
                view_data["details"] = [
                    {
                        "term_a": term_a.id,
                        "term_a_name": term_a.name,
                        "term_b": term_b.id,
                        "term_b_name": term_b.name,
                        "score": view_data["result"],
                    }
                ]

        elif mode == "gene" and gene_a_query and gene_b_query:
            gene_a = self.repo.find_gene(gene_a_query)
            gene_b = self.repo.find_gene(gene_b_query)

            if not gene_a or not gene_b:
                view_data["error"] = "One or both genes not found."
            else:
                calc = self.gene_calculators[strategy_name]
                view_data["result"], view_data["details"] = calc.calculate_similarity(
                    gene_a, gene_b
                )

        elif mode == "matrix" and gene_list_query:
            gene_symbols = [s.strip() for s in gene_list_query.split(",") if s.strip()]
            genes = []
            missing = []
            for sym in gene_symbols:
                g = self.repo.find_gene(sym)
                if g:
                    genes.append(g)
                else:
                    missing.append(sym)

            if len(genes) < 2:
                view_data["error"] = "Please provide at least two valid genes."
            else:
                calc = self.gene_calculators[strategy_name]
                labels, matrix = calc.calculate_matrix(genes)

                max_val = np.max(matrix) if np.max(matrix) > 0 else 1.0

                # Best Pair logic
                mask = np.ones(matrix.shape, dtype=bool)
                np.fill_diagonal(mask, 0)
                masked_matrix = matrix * mask
                i, j = np.unravel_index(masked_matrix.argmax(), masked_matrix.shape)

                normalized_matrix = (
                    matrix / max_val if max_val > 0 else np.zeros_like(matrix)
                )

                view_data["matrix_data"] = {
                    "labels": labels,
                    "matrix": matrix.tolist(),
                    "normalized_matrix": normalized_matrix.tolist(),
                    "max_val": float(max_val),
                    "best_pair": (labels[i], labels[j]),
                    "best_score": matrix[i, j],
                }
                if missing:
                    view_data["error"] = (
                        f"Some genes were not found: {', '.join(missing)}"
                    )

        return view_data
