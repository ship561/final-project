import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.regex.Pattern;

import net.sf.javaml.clustering.mcl.MarkovClustering;
import net.sf.javaml.clustering.mcl.SparseMatrix;
import net.sf.javaml.core.*;



public class predictMotif {

	ArrayList<String[]> fileIN(File f, ArrayList<String[]> motifs) {           //reads file in
		//StringBuffer contents = new StringBuffer();
		BufferedReader reader = null;
		String[] input = new String[350];
		ArrayList list = new ArrayList();
		try {
			reader = new BufferedReader(new FileReader(f));
			String text = null;

			while (reader.readLine() != null) {
				text = reader.readLine();
				if(Pattern.matches("^>.*", text) && list.size() > 0) {
					String[] m = new String[list.size()];
					list.toArray(m);
					motifs.add(m);
					list = new ArrayList();
				} else {
					list.add(text);
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifs;
	}
	
	double score(weighed_matrix P1, weighed_matrix F2) {
		double sumAln=0;
		double sum=0;
		double temp=0;
		double maxA=0;
		double maxPrf=0;
		for(int j=0; j<P1.Prf.get("A").length; j++) {
			sumAln=0;
			for(int si=0; si<F2.basefreq.get("A").length; si++) {
				for (Enumeration e = P1.Prf.elements() ; e.hasMoreElements() ;) {
					sumAln += F2.basefreq.get(e.nextElement())[si]*P1.Prf.get(e.nextElement())[j]; //sum(f2(b,s(i))*Prf1(b,i))
				}
			}
			temp = P1.I(j)*sumAln;
			if (temp > maxA)
				maxA = temp;
		}
		for(int i=0; i<P1.L; i++) {
			for (Enumeration e = P1.Prf.elements() ; e.hasMoreElements() ;) {
				temp = F2.Prf.get(e.nextElement())[i];	//max Prf(b,i)
				if (temp > maxPrf)
					maxPrf=temp;
			}
			sum += P1.I(i)*maxPrf;
		}
		return( maxA/(F2.motif.length*sum));
	}
	
	double sim(weighed_matrix M1, weighed_matrix M2) {
		double similarity;
		similarity = (score(M1,M2) + score(M2,M1))/2;
		return similarity;
	}
	Graph simgraph(weighed_matrix[] motifs) {
		int numMotif = motifs.length;
		Graph G = new Graph(numMotif);
		for(int i=0; i<numMotif-1; i++) {
			for (int j=1; j< numMotif; j++) {
				double similarity = sim(motifs[i], motifs[j]);
				if (similarity > 1) {					//Similarity needs to be greater than a cutoff beta
					G.adjacency_matrix[i][j] = 1;
					G.edgevalues[i][j] = similarity;
				}
			}
		}
		return G;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		MarkovClustering m = new MarkovClustering();
		SparseMatrix SM;
		
		predictMotif inputmotifs = new predictMotif();
		File file = new File("motif.txt");
		Graph simGraph;
		ArrayList<String[]> listMotifs = new ArrayList<String[]>();
		listMotifs = inputmotifs.fileIN(file, listMotifs);
		weighed_matrix[] matrixController = new weighed_matrix[listMotifs.size()];
		for(int i=0; i< listMotifs.size(); i++) {
			matrixController[i]= new weighed_matrix(listMotifs.get(i));
		}
		simGraph = new Graph(listMotifs.size());
		simGraph = inputmotifs.simgraph(matrixController);
		SM = new SparseMatrix(simGraph.adjacency_matrix);
		SM = m.run(SM, .1, .1, .1, 1);
	}

}
