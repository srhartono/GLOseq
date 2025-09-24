import pydot
import codecs
import re
import sys
import os
import graphviz
#dot = graphviz.Digraph(comment='The Round Table')

def read_with_bom_detection(filename):
    """Read file with automatic BOM detection and encoding handling."""
    encodings = ['utf-16', 'utf-16-le', 'utf-16-be', 'utf-8-sig', 'utf-8']
    
    for encoding in encodings:
        try:
            with open(filename, 'r', encoding=encoding) as f:
                content = f.read()
                content = content.lstrip('\ufeff')  # Remove any remaining BOM
                return content
        except (UnicodeDecodeError, UnicodeError):
            continue
    
    # Fallback: binary mode with manual BOM handling
    with open(filename, 'rb') as f:
        raw_data = f.read()
        if raw_data.startswith(codecs.BOM_UTF16_LE):
            return raw_data[2:].decode('utf-16-le')
        elif raw_data.startswith(codecs.BOM_UTF16_BE):
            return raw_data[2:].decode('utf-16-be')
        elif raw_data.startswith(codecs.BOM_UTF8):
            return raw_data[3:].decode('utf-8')
        else:
            return raw_data.decode('utf-8', errors='ignore')

def extract_dot_graph(dotFile):
    """Extract the DOT graph portion from a file that may contain other text."""

    content = None
    digraph_match = None

    with open(dotFile, 'r',encoding='utf-8') as f:
        content = f.read()
        digraph_match = re.search(r'digraph.+', content, re.MULTILINE)
        if digraph_match:
            print(f"Found 'digraph' at {dotFile}:{digraph_match.start()} to {digraph_match.end()} ({digraph_match.group()})\n")
        else:
            print(f'Error: digraph not found in {dotFile}! content=\n{content}\n')
            sys.exit(1)

    clean_dot = content[digraph_match.start():]
    return(clean_dot)

def check_graphviz():
    """Check if Graphviz is installed and available."""
    import subprocess
    try:
        subprocess.run(['dot', '-V'], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def main(dotFile, dot_exe_location="dot"):
    dot_exe_location_windows = "C:/Program Files/graphviz/bin/dot.exe"

    if not os.path.exists(dot_exe_location):
        print(f"original dot executable '{dot_exe_location}' not found! Trying {dot_exe_location_windows} (windows)!")
        dot_exe_location = dot_exe_location_windows
    if not os.path.exists(dot_exe_location):
        print(f"Error: dot executable '{dot_exe_location}' (windows) also not found! Stopping!")
        sys.exit(1)

    # dotFile = "dag.dot"
    output_file = "dag.pdf"
    
    if not os.path.exists(dotFile):
        print(f"Error: {dotFile} not found!")
        sys.exit(1)
    
    try:
        # Read and parse the DOT file
        print(f"Reading {dotFile}...")
        clean_dot = extract_dot_graph(dotFile)
        # with open(dotFile, 'r') as f:
        #graphs = pydot.graph_from_dot_file(clean_dotFile)
       
        graphs = pydot.graph_from_dot_data(clean_dot)
        if not graphs or len(graphs) == 0:
            print("Error: Failed to parse DOT file - no valid graph found")
            sys.exit(1)
            
        graph = graphs[0]
        print(f"Successfully parsed DOT file!")
        print(f"Graph has {len(graph.get_nodes())} nodes and {len(graph.get_edges())} edges")
        
        print(f"Generating {output_file}...")
        graph.write_pdf(output_file,prog=dot_exe_location)

        print(f"Successfully generated {output_file}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
