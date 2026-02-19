# Gene Ontology Analysis System
A web-based tool for analyzing Gene Ontology (GO) terms and gene annotations. This application allows researchers to explore the ontology structure, search for genes, and calculate semantic similarity scores.


## Installation

Install the required libraries directly using pip:

```
pip install flask pandas numpy
```
1. Open a terminal/cmd in the project folder.
2. Run the application:

```
python main.py
```

3. Once the server starts, open your web browser and navigate to:
   
   **http://127.0.0.1:5000**

   
## How to Use

**Upload Data**: The system requires data to function. Go to the **Upload** page and provide:
   - An OBO file (e.g., `go-basic.obo`)
   - A GAF file (e.g., `goa_human.gaf`)

Then the app will reload with the data, and various tools can be used. like:

- **Search**: Look up specific Genes (by symbol or name) or GO Terms (by ID or name).
- **Similarity**: Calculate similarity scores between two terms or two genes.
- **Pathfinder**: Find the path and distance between two GO terms.
- **Annotation Check**: Check if there is a gene annotation for a certain term.
