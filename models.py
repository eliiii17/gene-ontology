class GOTerm:
	"""
	Represents a single term in the Gene Ontology.
	Handles all the fields comprehensively.
	"""
	def __init__(self, id_, name, namespace):
		self.__id = id_
		self.__name = name[0].upper() + name[1:]
		self.__namespace = namespace.replace('_', ' ').title()

		self.__definition = ''
		self.__comments = []
		self.__alt_ids = []
		self.__synonyms = []
		self.__xrefs = []

		# Storing parent_id -> relationship_type (is_a, part_of, ...)
		# Using is_a and relationship fields in GO
		self.__parents = {}
		# The same as parents, only in reverse direction
		# This direction is not directly accessible, so it should be computed later
		self.__children = {}

	@property
	def id(self):
		return self.__id

	@property
	def name(self):
		return self.__name

	@property
	def namespace(self):
		return self.__namespace

	@property
	def definition(self):
		return self.__definition

	@property
	def comments(self):
		return self.__comments

	@property
	def alt_ids(self):
		return self.__alt_ids

	@property
	def synonyms(self):
		return self.__synonyms

	@property
	def xrefs(self):
		return self.__xrefs

	@property
	def parents(self):
		return self.__parents

	@property
	def children(self):
		return self.__children

	def set_def(self, definition):
		self.__definition = definition

	def add_comment(self, comment):
		self.__comments.append(comment)

	def add_alt_id(self, alt_id):
		self.alt_ids.append(alt_id)

	def add_synonym(self, synonym):
		self.__synonyms.append(synonym)

	def add_xref(self, xref):
		self.__xrefs.append(xref)

	def add_parent(self, parent_id, rel_type):
		self.__parents[parent_id] = rel_type

	def add_child(self, child_id, rel_type):
		self.__children[child_id] = rel_type

	def __repr__(self) -> str:
		return f"<GOTerm {self.id}: {self.name}>"


class Gene:
	"""
	Represents a human gene/protein.
	Aggregates multiple annotations from the GAF file into a single object.
	"""

	evidence_description_map = {
		"EXP": "Experiment",
		"IDA": "Direct Assay",
		"IPI": "Physical Interaction",
		"IMP": "Mutant Phenotype",
		"IGI": "Genetic Interaction",
		"IEP": "Expression Pattern",

		"HTP": "High Throughput Experiment",
		"HDA": "High Throughput Direct Assay",
		"HMP": "High Throughput Mutant Phenotype",
		"HGI": "High Throughput Genetic Interaction",
		"HEP": "High Throughput Expression Pattern",

		"IBA": "Biological aspect of Ancestor",
		"IBD": "Biological aspect of Descendant",
		"IKR": "Key Residues",
		"IRD": "Rapid Divergence",

		"IEA": "Electronic Annotation",

		"ISS": "Sequence or structural Similarity",
		"ISO": "Sequence Orthology",
		"ISA": "Sequence Alignment",
		"ISM": "Sequence Model",
		"IGC": "Genomic Context",
		"RCA": "Reviewed Computational Analysis",

		"TAS": "Traceable Author Statement",
		"NAS": "Non-traceable Author Statement",
		"IC": "Inferred by Curator",
		"ND": "No biological Data available"
	}

	def __init__(self, id_, symbol, name):
		# self._db = db  # column 1 in GAF
		self.__id = id_  # column 2
		self.__symbol = symbol  # column 3
		self.__name = name  # column 10

		self.__synonyms = []  # column 11

		# functional profile, the annotations
		# each annotaion is like this:
		# {'go_id': str, 'evidence': str, 'aspect': str, 'relation' str, 'count': int}
		self.__annotations = []

	def add_synonym(self, synonym):
		self.__synonyms.append(synonym)

	def add_annotation(self, annotation):
		"""
		Adds a new annotation or updates an existing one if the functional
		details (GO ID, Evidence, Aspect, Relation) are identical.
		Aggregates references and increments a count.
		"""
		# The fields that define uniqueness for an annotation 'type'
		key_fields = ['go_id', 'evidence', 'aspect', 'relation']

		found_ann = None
		for existing in self.annotations:
			# Check if all key fields match
			if all(existing[k] == annotation[k] for k in key_fields):
				found_ann = existing
				break

		if found_ann:
			found_ann['count'] += 1
		else:
			annotation['count'] = 1
			self.__annotations.append(annotation)

	@property
	def id(self):
		return self.__id

	@property
	def symbol(self):
		return self.__symbol

	@property
	def name(self):
		return self.__name

	@property
	def synonyms(self):
		return self.__synonyms

	@property
	def annotations(self):
		return self.__annotations

	@property
	def total_annotations_count(self):
		"""Returns the total number of annotation lines aggregated."""
		return sum(ann['count'] for ann in self.annotations)

	def get_valid_go_ids(self):
		"""Returns GO IDs, excluding the ones with 'NOT' relation."""
		return [
			annotation["go_id"]
			for annotation in self.annotations
			if "NOT" not in annotation["relation"]
		]

	def __repr__(self):
		return f"<Gene {self.id}: {self.symbol}: {len(self.__annotations)} annotations>"


	def __eq__(self, other):
		if not isinstance(other, Gene):
			raise NotImplementedError("cannot compare Gene with another type!")

		return self.id == other.id

	def __hash__(self):
		return hash(self.id)
