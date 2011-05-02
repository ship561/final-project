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

			while ((text=reader.readLine()) != null) {
				//text = reader.readLine();
				if(Pattern.matches("^>.*", text) && list.size() > 0) {
					String[] m = new String[list.size()];
					list.toArray(m);
					motifs.add(m);
					list = new ArrayList();
				} else if (! Pattern.matches("^>.*", text)){
					list.add(text);
				}
			}
			String[] m = new String[list.size()];
			list.toArray(m);
			motifs.add(m);
			reader.close();
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
		for(int j=0; j<P1.L; j++) {
			sumAln=0;
			for(int si=0; si<F2.L; si++) {
				for (Enumeration<String> e = P1.Prf.keys() ; e.hasMoreElements() ;) {
					String ele = e.nextElement();
					sumAln += F2.basefreq.get(ele)[si]*P1.Prf.get(ele)[j]; //sum(f2(b,s(i))*Prf1(b,i))
				}
			}
			temp = P1.I(j)*sumAln;
			if (temp > maxA)
				maxA = temp;
		}
		for(int i=0; i<P1.L; i++) {
			for (Enumeration<String> e = P1.Prf.keys() ; e.hasMoreElements() ;) {
				String b = e.nextElement();				//base b
				temp = F2.Prf.get(b)[i];	//max Prf(b,i)
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
		for(int i=0; i<numMotif; i++) {
			for (int j=0; j< numMotif; j++) {
				double similarity = sim(motifs[i], motifs[j]);
				if (similarity > .25) {					//Similarity needs to be greater than a cutoff beta
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
		for(int i=0; i< listMotifs.size(); i++) 
			matrixController[i]= new weighed_matrix(listMotifs.get(i));					//motifs are put into a controller array 
		simGraph = new Graph(listMotifs.size());										//creates a graph of the correct size
		simGraph = inputmotifs.simgraph(matrixController);								//adds edges
		SM = new SparseMatrix(simGraph.adjacency_matrix);
		SM = m.run(SM, .1, .1, .1, 1);
	}

}
