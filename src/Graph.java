
public class Graph {
	double[][] adjacency_matrix;
	double[][] edgevalues;
	
	Graph(int size) {
		adjacency_matrix 	= new double[size][size];
		edgevalues			= new double[size][size];
	}
	
	void addedge(int i, int j) {
		this.adjacency_matrix[i][j]=1;
		this.edgevalues[i][j]=1;
	}
}
