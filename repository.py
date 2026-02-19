import pandas as pd


class AnnotationRepository:
	"""
	Manages the collection of Gene objects and has methods
	for the UI and analysis.
	"""

	def __init__(self):
		# dict of Gene id: Gene object
		self.__genes = {}

		# dict, where each key is the GOTerm id, and each value
		# is a set of all the related Gene objects
		self.__term_gene_map = {}

	def _construct_term_gene_map(self):
		"""
		Creates a reverse mapping (GO Term -> Genes) for faster access.
		This allows instant retrieval of all genes for a specific GO term.
		"""
		for gene in self.__genes.values():
			# filtering NOT relations for accurate mapping
			go_ids = gene.get_valid_go_ids()

			for go_id in go_ids:
				if go_id not in self.__term_gene_map:
					self.__term_gene_map[go_id] = set()
				self.__term_gene_map[go_id].add(gene)

	def load_from_parser(self, parser):
		"""
		Loads genes by using the parser and saves genes.
		Also constructs the term to gene map(dict).
		parser is a parser object that has a .parse method
		"""
		self.__genes = parser.parse()

		self._construct_term_gene_map()

	def get_gene(self, gene_id):
		"""Gets a gene by its exact ID. Returns None if not found."""
		return self.__genes.get(gene_id)

	def find_gene(self, query):
		"""
		Finds a gene by ID or symbol.
		First attempts an exact ID match, then searches by symbol/name.
		Returns the first match found, or None.
		"""
		# Try exact ID match first
		gene = self.get_gene(query)
		if gene:
			return gene

		# If not found, search
		results = self.search_gene(query)
		if results:
			return results[0]

		return None

	@property
	def genes(self):
		return list(self.__genes.values())

	def search_gene(self, search_str, limit=None):
		"""
		Search by symbol, synonym, or name with prioritization.
		Priority: 1. Exact Symbol, 2. Exact Synonym, 3. Partial Symbol,
		4. Partial Synonym, 5. Partial Name.
		"""
		search_str = search_str.strip().lower()
		if not search_str:
			return []
		if limit and limit > len(self.__genes):
			limit = None

		exact_symbol_matches = []
		exact_synonym_matches = []
		partial_symbol_matches = []
		partial_synonym_matches = []
		partial_name_matches = []

		for gene in self.genes:
			symbol_lower = gene.symbol.lower()
			name_lower = gene.name.lower()
			synonyms_lower = [synonym.lower() for synonym in gene.synonyms]

			if search_str == symbol_lower:
				exact_symbol_matches.append(gene)
			elif search_str in synonyms_lower:
				exact_synonym_matches.append(gene)
			elif search_str in symbol_lower:
				partial_symbol_matches.append(gene)
			elif any(search_str in synonym for synonym in synonyms_lower):
				partial_synonym_matches.append(gene)
			elif search_str in name_lower:
				partial_name_matches.append(gene)

		results =  exact_symbol_matches + exact_synonym_matches + partial_symbol_matches + partial_synonym_matches + partial_name_matches

		return results[:limit]

	def to_dataframe(self):
		"""	Exports all annotations to a Pandas DataFrame."""
		data = []
		for gene in self.__genes.values():
			for ann in gene.annotations:
				data.append({
					"gene_id": gene.id,
					"go_id": ann["go_id"],
					"evidence": ann["evidence"],
					"aspect": ann["aspect"],
					"relation": ann["relation"],
					"count": ann["count"]
				})
		return pd.DataFrame(data)

	def get_annotated_terms(self):
		"""Returns a list of all GO IDs that have at least one gene annotation."""
		return list(self.__term_gene_map.keys())

	def get_genes_for_term(self, go_id):
		"""Returns all genes explicitly annotated to a specific GO Term as a list."""
		return list(self.__term_gene_map.get(go_id, set()))

	def get_genes_for_term_recursive(self, go_id, graph, relations=["is_a"]):
		"""Returns genes annotated to the specific term and all its children as a list."""

		results_set = set(self.get_genes_for_term(go_id))

		children_ids_set = graph.get_descendants(go_id, relations=relations)

		for child_id in children_ids_set:
			if child_id in self.__term_gene_map:
				results_set.update(self.__term_gene_map[child_id])

		results = list(results_set)
		results.sort(key=lambda gene: gene.symbol)  # to make results consistent

		return results

	def check_gene_term_annotation(self, gene_query, term_id, graph):
		"""
		Checks if a gene is annotated to a term or any of its descendants.
		Returns a dict like:
		{
			"is_annotated": True/False,
			"gene": Gene,
			"term": GOTerm,
			"other_genes": list of other genes,
			"error": str		
		}
		"""

		gene = self.find_gene(gene_query)
		if not gene:
			return {"error": f"Gene '{gene_query}' not found."}

		term = graph.get_term(term_id)
		if not term:
			return {"error": f"Term '{term_id}' not found."}

		# all descendants of the term
		term_family = graph.get_descendants(term_id, relations=None)

		gene_annotations = set(gene.get_valid_go_ids())
		is_annotated = len(gene_annotations.intersection(term_family)) > 0

		if not is_annotated:
			return {"gene": gene, "term": term, "is_annotated": False}
		else:
			# if annotated, find all other genes in that term family.
			all_genes = self.get_genes_for_term_recursive(
				term_id, graph, relations=None
			)
			other_genes = [g for g in all_genes if g.id != gene.id]

			return {
				"gene": gene,
				"term": term,
				"is_annotated": True,
				"other_genes": other_genes
			}
