import os

from flask import Flask, flash, redirect, render_template, request, url_for

from analysis import InformationContentCalculator, SimilarityService
from statistics import StatisticsAnalyzer
from models import Gene
from ontology import OntologyGraph
from parsers import GAFParser, OBOParser
from repository import AnnotationRepository

app = Flask("Gene Ontology")
app.secret_key = "1"

# Global variables for data
graph = None
repo = None
ic_calc = None
similarity_service = None
system_ready = False


def load_data(obo_path="data/go-basic.obo", gaf_path="data/goa_human.gaf"):
    """Loads or reloads the ontology and annotation data."""
    global graph, repo, ic_calc, strategies, gene_sim_calculators, system_ready, similarity_service

    print(f"Loading ontology from {obo_path}... ", end="", flush=True)
    if os.path.exists(obo_path):
        graph = OntologyGraph()
        graph.load_from_parser(OBOParser(obo_path))
        print("Done.")
    else:
        print("File not found.")

    print(f"Loading annotations from {gaf_path}... ", end="", flush=True)
    if os.path.exists(gaf_path):
        repo = AnnotationRepository()
        repo.load_from_parser(GAFParser(gaf_path))
        print("Done.")
    else:
        print("File not found.")

    system_ready = bool(repo and graph)

    if system_ready:
        print("Initializing analysis tools... ", end="", flush=True)
        ic_calc = InformationContentCalculator(graph, repo)
        similarity_service = SimilarityService(graph, repo, ic_calc)
        print("Done.")

        print(
            f"System Ready. Loaded {len(repo.genes)} genes and {len(graph.all_term_ids())} terms."
        )

# Initial load
load_data()


@app.route("/")
def index():
    """Dashboard Homepage"""

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    stats = {
        "gene_count": len(repo.genes),
        "term_count": len(graph.all_term_ids()),
    }

    return render_template("index.html", stats=stats, active_page="index")


@app.route("/upload", methods=["GET", "POST"])
def upload():
    """Handles file uploads for OBO and GAF files."""
    if request.method == "POST":
        os.makedirs("data", exist_ok=True)
        
        obo_file = request.files.get("obo_file")
        gaf_file = request.files.get("gaf_file")

        updated = False
        if obo_file and obo_file.filename:
            obo_path = os.path.join("data", "go-basic.obo")
            obo_file.save(obo_path)
            updated = True

        if gaf_file and gaf_file.filename:
            gaf_path = os.path.join("data", "goa_human.gaf")
            gaf_file.save(gaf_path)
            updated = True

        if updated:
            load_data()
            flash("Data files updated and system reloaded successfully!", "success")
            return redirect(url_for("index"))
        else:
            flash("No files were uploaded.", "warning")
            return render_template("upload.html", active_page="upload")

    return render_template("upload.html", active_page="upload")


@app.route("/stats")
def stats():
    """Displays statistical summaries of the loaded data."""

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    analyzer = StatisticsAnalyzer(graph, repo)
    stats_data = analyzer.calculate_all()

    return render_template("stats.html", stats=stats_data, active_page="stats")


@app.route("/search")
def search():
    """Renders the search page and handles search queries."""

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    query = request.args.get("query")
    search_type = request.args.get("search_type", "gene")
    results = None

    if query:
        if search_type == "gene":
            results = repo.search_gene(query)
        else:
            results = graph.search_term(query)

    return render_template(
        "search.html",
        results=results,
        query=query,
        search_type=search_type,
        active_page="search",
    )


@app.route("/similarity")
def similarity():
    """Handles the similarity calculation."""

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    view_data = similarity_service.process_request(request.args)
    view_data["active_page"] = "similarity"

    return render_template("similarity.html", **view_data)


@app.route("/pathfinder")
def pathfinder():
    """Renders the pathfinder page and shows direct or common ancestor paths."""

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    term_a_id = request.args.get("term_a_id")
    term_b_id = request.args.get("term_b_id")

    # These are kept to repopulate the form after submission
    term_a_query = request.args.get("term_a_query")
    term_b_query = request.args.get("term_b_query")

    path_data = {
        "term_a_query": term_a_query,
        "term_b_query": term_b_query,
        "term_a_id": term_a_id,
        "term_b_id": term_b_id,
        "path_type": None,
        "active_page": "pathfinder",
    }

    if term_a_id and term_b_id:
        result = graph.get_relationship_between_terms(term_a_id, term_b_id)
        if result:
            path_data.update(result)
        else:
            path_data['error'] = "No direct path or common ancestor found."
    return render_template("pathfinder.html", **path_data)


@app.route("/annotation-check")
def annotation_check():
    """
    Handles the gene annotation check. The user provides a gene and a term,
    and the system checks if the gene is annotated to that term or any of its
    children. Autocomplete is used on the frontend, so this route handles the
    final calculation.
    """

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    gene_query = request.args.get("gene_query")
    term_query = request.args.get("term_query")  # Kept for repopulating form
    term_id = request.args.get("term_id")

    result_data = None
    if gene_query and term_id:
        result_data = repo.check_gene_term_annotation(gene_query, term_id, graph)

    return render_template(
        "annotation_check.html",
        gene_query=gene_query,
        term_query=term_query,
        term_id=term_id,
        result=result_data,
        active_page="annotation_check",
    )


@app.route("/term/<term_id>")
def term_detail(term_id):
    """Displays the details for a specific GO term."""

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    term = graph.get_term(term_id)
    if not term:
        flash(f"Term {term_id} not found!", "error")
        return redirect(request.referrer or url_for("index"))

    neighborhood = graph.get_neighborhood(term_id, relations=None)
    genes = repo.get_genes_for_term_recursive(term_id, graph, relations=None)

    return render_template(
        "term_detail.html",
        term=term,
        neighborhood=neighborhood,
        genes=genes,
        graph=graph,
    )


@app.route("/gene/<gene_id>")
def gene_detail(gene_id):
    """Displays the details for a specific gene."""

    if not system_ready:
        flash("Data not loaded. Please upload files.", "error")
        return redirect(url_for("upload"))

    gene = repo.get_gene(gene_id)
    if not gene:
        flash(f"Gene {gene_id} not found!", "error")
        return redirect(request.referrer or url_for("index"))

    # Get full term objects for annotations to display names
    term_names = {}

    # Calculate summaries
    aspect_counts = {"P": 0, "F": 0, "C": 0}

    for an in gene.annotations:
        term = graph.get_term(an["go_id"])
        if term:
            term_names[an["go_id"]] = term.name

        # Summarize aspects
        if an["aspect"] in aspect_counts:
            aspect_counts[an["aspect"]] += an["count"]

    # Sort annotations by count to find top terms
    sorted_anns = sorted(gene.annotations, key=lambda x: x["count"], reverse=True)
    top_annotations = sorted_anns[:5]

    return render_template(
        "gene_detail.html",
        gene=gene,
        annotations=gene.annotations,
        term_names=term_names,
        evidence_description_map=Gene.evidence_description_map,
        aspect_counts=aspect_counts,
        top_annotations=top_annotations,
    )


@app.route("/api/search-terms")
def api_search_terms():
    """API endpoint to search for GO terms and return them as JSON."""

    if not system_ready:
        return {}

    query = request.args.get("q", "")
    if len(query) < 3:
        # Don't search for very short strings to avoid too many results
        return {"results": []}

    terms = graph.search_term(query)

    # Limit the number of results to keep the dropdown manageable
    results = [
        {
            "id": term.id,
            "name": term.name,
            "namespace": term.namespace.split()[1].upper()[0],  # only F, P, C
        }
        for term in terms[:15]  # Return at most 15 results
    ]
    return {"results": results}


@app.route("/api/search-genes")
def api_search_genes():
    """API endpoint to search for genes and return them as JSON."""

    if not system_ready:
        return {}

    query = request.args.get("q", "")
    if len(query) < 2:
        return {"results": []}

    genes = repo.search_gene(query, limit=15)
    results = [
        {"id": gene.id, "symbol": gene.symbol, "name": gene.name} for gene in genes
    ]
    return {"results": results}


@app.route("/api/top-genes")
def api_top_genes():
    """Returns the top 10 annotated genes as a comma-separated string."""
    if not repo:
        return {"genes": ""}

    # Calculate top genes dynamically
    anns_df = repo.to_dataframe()
    top_genes = (
        anns_df.groupby("gene_id")["count"].sum().sort_values(ascending=False).head(10)
    )

    symbols = []
    for gene_id in top_genes.index:
        gene = repo.get_gene(gene_id)
        if gene:
            symbols.append(gene.symbol)

    return {"genes": ", ".join(symbols)}


app.run()
