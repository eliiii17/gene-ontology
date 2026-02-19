class StatisticsAnalyzer:
    """
    Encapsulates logic for calculating statistical summaries of the
    Ontology and Annotation data.
    """

    def __init__(self, graph, repo):
        self.graph = graph
        self.repo = repo
        self.graph_df = graph.to_dataframe()
        self.anns_df = repo.to_dataframe()

    def calculate_all(self):
        """Returns a dictionary containing all statistical metrics."""
        return {
            **self._ontology_structure_stats(),
            **self._annotation_general_stats(),
            **self._top_rankings(),
            **self._specificity_metrics(),
            **self._aspect_stats(),
            **self._evidence_stats(),
        }

    def _ontology_structure_stats(self):
        return {
            "ns_counts": self.graph_df["namespace"].value_counts().to_dict(),
            "avg_parents": round(self.graph_df["parents_count"].mean(), 2),
            "avg_children": round(self.graph_df["children_count"].mean(), 2),
        }

    def _annotation_general_stats(self):
        anns_per_gene = self.anns_df.groupby("gene_id")["count"].sum()
        return {
            "avg_anns_per_gene": round(anns_per_gene.mean(), 2),
            "total_annotations": int(self.anns_df["count"].sum())
        }

    def _top_rankings(self):
        # Top 10 genes
        anns_per_gene = self.anns_df.groupby("gene_id")["count"].sum()
        top_genes = anns_per_gene.sort_values(ascending=False).head(10)

        top_genes_data = []
        for gene_id, count in top_genes.items():
            gene = self.repo.get_gene(gene_id)
            if gene:
                top_genes_data.append(
                    {"id": gene.id, "symbol": gene.symbol, "count": count}
                )

        # Top 10 terms
        term_counts = self.anns_df.groupby("go_id")["count"].sum()
        top_terms = term_counts.sort_values(ascending=False).head(10)

        top_terms_data = []
        for term_id, count in top_terms.items():
            term = self.graph.get_term(term_id)
            name = term.name if term else "Unknown"
            top_terms_data.append({"id": term_id, "name": name, "count": int(count)})

        return {"top_genes": top_genes_data, "top_terms": top_terms_data}

    def _specificity_metrics(self):
        total_depth_sum = 0
        total_annotations_count = 0
        leaf_annotations_count = 0

        # Weighted calculation based on term usage count
        term_counts = self.anns_df.groupby("go_id")["count"].sum()

        for term_id, count in term_counts.items():
            term = self.graph.get_term(term_id)
            if term:
                d = self.graph.get_depth(term_id)
                total_depth_sum += d * count
                total_annotations_count += count

                if not term.children:
                    leaf_annotations_count += count

        avg_depth = (
            total_depth_sum / total_annotations_count if total_annotations_count else 0
        )
        leaf_ratio = (
            (leaf_annotations_count / total_annotations_count * 100)
            if total_annotations_count
            else 0
        )

        return {
            "avg_annotation_depth": round(avg_depth, 2),
            "leaf_annotation_ratio": round(leaf_ratio, 1),
        }

    def _aspect_stats(self):
        aspect_counts = self.anns_df.groupby("aspect")["count"].sum().to_dict()
        aspect_map = {
            "P": "Biological Process",
            "F": "Molecular Function",
            "C": "Cellular Component",
        }
        return {
            "aspect_stats": {
                aspect_map[k]: int(v) for k, v in aspect_counts.items()
            }
        }

    def _evidence_stats(self):
        evidence_counts = (
            self.anns_df.groupby("evidence")["count"]
            .sum()
            .sort_values(ascending=False)
            .head(10)
            .to_dict()
        )
        return {"evidence_counts": evidence_counts}
