from collections import deque
import pandas as pd

class OntologyGraph:
	"""
	Holds all the terms of the gene ontology and manages the graph structure
	based on different relations like is_a, part_of, ...
	"""

	def __init__(self):
		# for storing all the terms. dict that maps GO ids to GOTerm objects
		self.__terms = {}
		self.__depth_cache = {}

	def get_term(self, term_id):
		if term_id not in self.__terms:
			print("Invalid ID! Does not exist")
			return None

		return self.__terms[term_id]

	def all_term_ids(self):
		# using list to keep order
		return list(self.__terms.keys())

	def load_from_parser(self, parser):
		"""
		Loads terms by using the parser and reconstructs graph edges.
		parser is an object with a .parse method
		"""

		self.__terms = parser.parse()  # this will only create parents for each term
		# and children is still empty

		self._calculate_children()

	def _calculate_children(self):
		"""
		Iterates over all terms to infer and calculate children attributes.
		"""
		for term_id, term in self.__terms.items():
			for parent_id, rel_type in term.parents.items():
				if parent_id in self.__terms:
					self.__terms[parent_id].add_child(term_id, rel_type)

	def find_downward_path(self, start_id, end_id):
		"""
		Finds the shortest path from an ancestor (start_id) to a descendant (end_id).
		This search only goes down the graph (via children).
		Returns a list of dicts like {'term': GOTerm, 'relation': str} in order
		Returns None if no path
		"""
		if start_id not in self.__terms or end_id not in self.__terms:
			return None

		queue = deque([[start_id]])
		visited = {start_id}

		while queue:
			path = queue.popleft()
			current_id = path[-1]

			if current_id == end_id:
				# Path found
				enriched_path = []
				for i in range(len(path)):
					node_id = path[i]
					term = self.get_term(node_id)
					relation_info = {"term": term}

					if i + 1 < len(path):
						next_node_id = path[i + 1]
						relation_info["relation"] = term.children.get(next_node_id)

					enriched_path.append(relation_info)
				return enriched_path

			term = self.get_term(current_id)
			for child_id in term.children:
				if child_id not in visited:
					visited.add(child_id)
					new_path = list(path)
					new_path.append(child_id)
					queue.append(new_path)

		return None  # No downward path found

	def find_lca_path(self, term_a_id, term_b_id):
		"""
		Finds the Lowest(nearest) Common Ancestor (LCA) and the paths to it.
		Returns a dict with the lca, and the two paths from it.
		like 
		{
			"lca": GOTerm,
			"path_a": the path output of self.find_downward_path method, 
			"path_b": the path output of self.find_downward_path method 
		}
		returns None if not found
		"""
		ancestors_a = self.get_ancestors(term_a_id)
		ancestors_b = self.get_ancestors(term_b_id)

		common_ancestors = ancestors_a.intersection(ancestors_b)
		if not common_ancestors:
			return None

		# Find the LCA (nearest common ancestor)
		lca = max(common_ancestors, key=lambda term_id: self.get_depth(term_id))

		# Get the paths from the LCA down to each term
		path_a = self.find_downward_path(lca, term_a_id)
		path_b = self.find_downward_path(lca, term_b_id)

		if not path_a or not path_b:
			return None

		return {"lca": self.get_term(lca), "path_a": path_a, "path_b": path_b}

	def get_relationship_between_terms(self, term_a_id, term_b_id):
		"""
		Determines the relationship between two terms, finding either a direct
		ancestral path or a path via the Lowest Common Ancestor (LCA).
		Returns a dictionary structured for the template(path finder).
		{
			"path_type": 'direct' or 'lca',
			'path': path data,
			'lca_path': lca path data,
			"path_direction": str
		}
		returns None if no path found.
		"""
		term_a = self.get_term(term_a_id)
		term_b = self.get_term(term_b_id)

		if not term_a or not term_b:
			return {"error": "One or both term IDs are invalid."}

		# Try direct path A -> B
		direct_path = self.find_downward_path(term_a_id, term_b_id)
		if direct_path:
			return {
				"path_type": "direct",
				"path": direct_path,
				"path_direction": f"Path from {term_a.name} to {term_b.name}",
			}

		# Try direct path B -> A
		direct_path = self.find_downward_path(term_b_id, term_a_id)
		if direct_path:
			return {
				"path_type": "direct",
				"path": direct_path,
				"path_direction": f"Path from {term_b.name} to {term_a.name}",
			}

		# If no direct path, find common ancestor
		lca_data = self.find_lca_path(term_a_id, term_b_id)
		if lca_data:
			return {"path_type": "lca", "lca_path": lca_data}

		return None

	def search_term(self, query):
		"""
		Searches for terms by ID or name, prioritizing exact ID matches, then
		partial ID matches, then partial name matches.
		"""
		query = query.strip().lower()

		exact_id_matches = []
		exact_name_matches = []
		exact_synonym_matches = []
		partial_id_matches = []
		partial_name_matches = []
		partial_synonym_matches = []

		for term in self.__terms.values():
			id_lower = term.id.lower()
			name_lower = term.name.lower()
			synonyms_lower = [syn.lower() for syn in term.synonyms]

			if query == id_lower:
				exact_id_matches.append(term)
			elif query in id_lower:
				partial_id_matches.append(term)
			elif query == name_lower:
				exact_name_matches.append(term)
			elif query in synonyms_lower:
				exact_synonym_matches.append(term)
			elif query in name_lower:
				partial_name_matches.append(term)
			elif any(query in synonym for synonym in synonyms_lower):
				partial_synonym_matches.append(term)

		return exact_id_matches + partial_id_matches + exact_name_matches + exact_synonym_matches + partial_name_matches + partial_synonym_matches

	def get_ancestors(self, term_id, relations=["is_a"]):
		"""
		relations: is a list of relation types to follow, like is_a, part_of, ...
		if it is None, we follow every relation
		Returns a set of ancestor IDs.
		"""

		if term_id not in self.__terms:
			return set()

		ancestors = set()

		# BFS algorithm to find all ancestors
		queue = [term_id]
		while len(queue) > 0:
			current_id = queue.pop(0)

			if current_id not in ancestors:
				ancestors.add(current_id)

				if current_id in self.__terms:
					term = self.__terms[current_id]
				else:
					continue
				# checking parents to decide which edges to follow
				for parent_id, rel_type in term.parents.items():
					if (relations is None) or (rel_type in relations):
						queue.append(parent_id)

		return ancestors

	def get_descendants(self, term_id, relations=["is_a"]):
		"""
		relations: is a list of relation types to follow, like is_a, part_of, ...
		if it is None, we follow every relation
		Returns a set of descendant IDs.
		"""
		# BFS algorithm to find all descendants
		if term_id not in self.__terms:
			return set()

		descendants = set()

		queue = [term_id]

		while len(queue) > 0:
			current_id = queue.pop(0)

			if current_id not in descendants:
				descendants.add(current_id)

				if current_id in self.__terms:
					term = self.__terms[current_id]
				else:
					continue

				for child_id, rel_type in term.children.items():
					if (relations is None) or (rel_type in relations):
						queue.append(child_id)

		return descendants

	def get_depth(self, term_id):
		"""
		Calculates the depth of a term.
		The depth is the longest path from a root node.
		Uses cache for performance.
		"""
		if term_id in self.__depth_cache:
			return self.__depth_cache[term_id]

		term = self.get_term(term_id)
		if not term:
			return 0

		if not term.parents:
			self.__depth_cache[term_id] = 0
			return 0

		parent_depths = [self.get_depth(parent_id) for parent_id in term.parents]
		depth = max(parent_depths) + 1
		self.__depth_cache[term_id] = depth
		return depth

	def get_neighborhood(self, term_id, relations=["is_a"]):
		"""
		Returns immediate parents and children of a term.
		relations: is a list of relation types to follow, like is_a, part_of, ...
		if it is None, we follow every relation
		returns {"parents": parents dict, "children": children dict}
		"""
		term = self.get_term(term_id)
		if not term:
			print(f"{term_id} does not exist!")
			return {"parents": {}, "children": {}}

		if relations is None:
			return {"parents": term.parents, "children": term.children}

		# filter out only the desired relations
		parents = {
			p_id: r_type
			for p_id, r_type in term.parents.items()
			if r_type in relations
		}
		children = {
			c_id: r_type
			for c_id, r_type in term.children.items()
			if r_type in relations
		}

		return {"parents": parents, "children": children}

	def to_dataframe(self):
		"""
		Exports the ontology terms to a Pandas DataFrame.
		"""
		data = []
		for term in self.__terms.values():
			data.append(
				{
					"id": term.id,
					"name": term.name,
					"namespace": term.namespace,
					"definition": term.definition,
					"parents_count": len(term.parents),
					"children_count": len(term.children),
				}
			)
		return pd.DataFrame(data)
