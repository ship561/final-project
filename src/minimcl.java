import net.sf.javaml.clustering.mcl.MarkovClustering;
import net.sf.javaml.clustering.mcl.SparseMatrix;


public class minimcl {
	double[][] addloops(double[][] adjacency_matrix) {
		int n= adjacency_matrix.length;
		for(int i=0; i<n; i++ ) {
			double max=0;
			for(int j=0; j<n; j++) { 
				for(int k=0; k<n; k++) {
					if(adjacency_matrix[i][k] > max) 
						max = adjacency_matrix[i][k];
				}
				if(max == 0) 
					max=1;
				if(i==j) 
					adjacency_matrix[i][j] = max;
			}
		}
		return adjacency_matrix;
	}
	
	double[][] make_stochastic(double[][] matrix) {
		SparseMatrix SM = new SparseMatrix(matrix);
		SM.normaliseRows();
		matrix = SM.getDense();
		return matrix;
	}
	
	SparseMatrix mcl(SparseMatrix mx, double I) {
		double chaos=1;
		int i = 1;
		MarkovClustering MC = new MarkovClustering();
		while (chaos > .00001) {
			SparseMatrix sq;
			sq = MC.expand(mx);
			chaos = MC.inflate(sq, I, .0001);
			mx = sq;
			i++;
		}
		return mx;
	}
	
	String[] getClusters(SparseMatrix mx) {
		mx.prune(.01);
		mx.normaliseRows();
		int size = mx.getSize()[0];
		String[] str = new String[size];
		for(int col=0; col<size; col++) {
			str[col] = "";
			for(int row=0; row<size; row++) {
				if(mx.get(row,col) != 0) {
					str[col] += row + " ";
				}
			}
		}
		return str;
	}
}
