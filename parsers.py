import os
from abc import ABC, abstractmethod

import pandas as pd

from models import Gene, GOTerm


class BaseParser(ABC):
    """
    Abstract base class for all data parsers in the GO analysis system.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"{filepath} does not exist!")

    @abstractmethod
    def parse(self):
        pass


class OBOParser(BaseParser):
    """Parses Gene Ontology OBO files into GOTerm objects."""

    def parse(self):
        """
        Reads the OBO file and returns a dictionary of {GO_ID: GOTerm}.
        """

        # final dict to return
        terms = {}

        with open(self.filepath, "r") as file:
            term_lines = []  # each item is a line(str)
            inside_term_block = False

            for line in file:
                line = line.strip()
                # skip empty lines
                if not line:
                    continue

                if line.startswith("[Term]"):
                    if inside_term_block and term_lines:
                        # we have reached to the end of a Term block
                        term = self._process_term(term_lines)
                        if term:
                            terms[term.id] = term

                    # reset for next term
                    term_lines = []
                    inside_term_block = True
                elif line.startswith("[Typedef]"):
                    if inside_term_block and term_lines:
                        term = self._process_term(term_lines)
                        if term:
                            terms[term.id] = term

                    # reset for next term
                    term_lines = []
                    inside_term_block = False

                elif inside_term_block:
                    term_lines.append(line)

        return terms

    def _process_term(self, term_lines):
        """
        Processes a term block with term_lines which is a list of lines(str)
        returns the GOTerm object.
        returns None if invalid or obsolete
        """
        # setting default values in case of non-existing fields
        def_ = ""
        alt_ids = []
        synonyms = []
        parents = {}  # dict, term_id: rel_type
        xrefs = []
        comments = []
        namespace = "gene_ontology"  # default namesoace according to OBO file

        id_ = ""
        name = ""

        for line in term_lines:
            line = line.strip()

            if line.startswith("is_obsolete:"):
                return None
            elif line.startswith("id:"):
                id_ = line.split()[1]
            elif line.startswith("alt_id:"):
                alt_ids.append(line.split()[1])
            elif line.startswith("name:"):
                name = line.split(":")[1].strip()
            elif line.startswith("def:"):
                def_ = line.split('"')[1].strip()
            elif line.startswith("comment:"):
                comments.append(line.split(":")[1].strip())
            elif line.startswith("synonym:"):
                synonyms.append(line.split('"')[1].strip())
            elif line.startswith("namespace:"):
                namespace = line.split(":")[1].strip()
            elif line.startswith("is_a:"):
                # each parents dict entry is like parent_id: rel_type
                parents[line.split()[1]] = "is_a"
            elif line.startswith("relationship:"):
                rel_fields = line.removeprefix("relationship:").strip().split()
                parents[rel_fields[1]] = rel_fields[0]
            elif line.startswith("xref:"):
                xrefs.append(line.removeprefix("xref: ").strip())

        if not (id_ and name):
            # is not a valid GO term! returning None
            return None

        # creating the GOTerm object, then setting other information
        term = GOTerm(id_=id_, name=name, namespace=namespace)
        term.set_def(def_)

        for alt_id in alt_ids:
            term.add_alt_id(alt_id)

        for syn in synonyms:
            term.add_synonym(syn)

        for comment in comments:
            term.add_comment(comment)

        for xref in xrefs:
            term.add_xref(xref)

        for parent_id, rel_type in parents.items():
            term.add_parent(parent_id, rel_type)

        return term


class GAFParser(BaseParser):
    """
    Parses Gene Association Format (GAF) file.
    Uses Pandas for file reading, then converts data into the Gene model.
    """

    COL_NAMES = [
        "db",  # 0 - skip
        "id",  # 1 (ID)
        "symbol",  # 2 (symbol)
        "relation",  # 3 (NOT, ...)
        "GO_id",  # 4  (GOTerm ID)
        "reference",  # 5  - skip
        "evidence_code",  # 6  (EXP, IDA, ...)
        "with_from",  # 7  - skip
        "aspect",  # 8  (P, F, C) same as the namespace of obo file
        "name",  # 9  (Full Name)
        "synonym",  # 10 (Synonyms separated by |)
    ]

    def parse(self):
        """
        Reads the GAF file and returns a dictionary of Gene objects.
        Each entry in dict is gene id(str): Gene Object
        """

        genes = {}

        df = pd.read_csv(
            self.filepath,
            sep="\t",  # since GAF is tabular
            comment="!",  # skip header lines starting with !
            header=None,
            names=self.COL_NAMES,
            usecols=range(11),  # only the first 11 columns
            dtype=str
        )

        for row in df.itertuples(index=False):
            # create a new Gene if it is the first time in the data
            if row.id not in genes:
                name = row.name if pd.notna(row.name) else ""
                new_gene = Gene(id_=row.id, symbol=row.symbol, name=name)

                if pd.notna(row.synonym):
                    synonyms = str(row.synonym).split("|")
                    for s in synonyms:
                        new_gene.add_synonym(s)

                # add new gene to dict
                genes[row.id] = new_gene

            # creating the annotation record
            relation = row.relation if pd.notna(row.relation) else ""

            annotation = {
                "go_id": row.GO_id,
                "evidence": row.evidence_code,
                "aspect": row.aspect,
                "relation": relation
            }

            # add the annotaion
            genes[row.id].add_annotation(annotation)

        return genes
