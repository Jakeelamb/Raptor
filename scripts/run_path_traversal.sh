#!/bin/bash
set -e

# Path Traversal Example Script
echo "GFA Path Traversal Example Script"
echo "================================="

# Create directory for test data if it doesn't exist
mkdir -p test_data
cd test_data

# Create a sample GFA file with paths
cat > sample_graph.gfa << 'EOF'
H	VN:Z:1.0
S	1	ACGTACGTACGT
S	2	TTTTGGGGCCCC
S	3	AAAATTTTGGGG
S	4	CCCCAAAATTTT
P	transcript_1	1+,2+,3+	*
P	transcript_2	1+,3-,4+	*
P	transcript_3	2+,3+,4-	*
EOF

echo "Created sample GFA file with path definitions."

# Create segment sequences file
cat > segment_sequences.tsv << 'EOF'
1	ACGTACGTACGT
2	TTTTGGGGCCCC
3	AAAATTTTGGGG
4	CCCCAAAATTTT
EOF

echo "Created segment sequences file."

# Run path traversal
echo "Running path traversal..."
cd ..
cargo run --release -- traverse \
    -i test_data/sample_graph.gfa \
    -s test_data/segment_sequences.tsv \
    -o test_data/output \
    --formats fasta,dot,json,tsv \
    --visualize \
    --metadata

echo ""
echo "Generated output files:"
ls -l test_data/output*

echo ""
echo "Example FASTA content:"
head -n 5 test_data/output.fasta

echo ""
echo "Example metadata content:"
head -n 5 test_data/output.path_stats.tsv

echo ""
echo "To visualize the DOT graph, you can use Graphviz:"
echo "dot -Tsvg test_data/output.dot -o test_data/output.svg"

# Generate visualization if dot is installed
if command -v dot &> /dev/null; then
    echo "Graphviz detected, generating SVG visualization..."
    dot -Tsvg test_data/output.dot -o test_data/output.svg
    echo "Visualization saved to test_data/output.svg"
else
    echo "Graphviz not found. Install it to visualize the graph."
fi 