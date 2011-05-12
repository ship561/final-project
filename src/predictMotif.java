import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Pattern;

import net.sf.javaml.clustering.mcl.SparseMatrix;



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
		for(int si=0; si<F2.L; si++) {
			sumAln=0;
			int i=0;
			while(si+i < F2.L && i < P1.L) {
				for (Enumeration<String> e = P1.Prf.keys() ; e.hasMoreElements() ;) {
					String ele = e.nextElement();
					sumAln += F2.basefreq.get(ele)[si+i]*P1.Prf.get(ele)[i]*P1.I(i); //I(i,P1)*sum(f2(b,s(i))*Prf1(b,i))
				}
				i++;
			}
			//temp = P1.I(si)*sumAln;
			if (sumAln > maxA)
				maxA = sumAln;
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
		double beta = .05; 				//minimum value for making an edge on the graph combining 2 nodes
		Graph G = new Graph(numMotif);
		for(int i=0; i<numMotif; i++) {
			for (int j=0; j< numMotif; j++) {
				double similarity = sim(motifs[i], motifs[j]);
				if (similarity > beta) {					//Similarity needs to be greater than a cutoff beta
					G.adjacency_matrix[i][j] = 1;
					G.edgevalues[i][j] = similarity;
				}
			}
		}
		return G;
	}
	String checkclique(String node, List<String> c, Graph G) {
		List<String> clique = c;
		int n = clique.size();
		double numEdge=n+1, weight=0;
		while (n >= 3) {
			double[] min = {n+1,0,0}; //number edges, summed weight of edges, index to be deleted 
			int edge = 0;
			for(int i=0; i<n; i++) {
				int v1 = Integer.parseInt(clique.get(i));		//first vertex
				numEdge=0;
				weight=0;
				for(int j=0; j<n; j++) {
					int v2 = Integer.parseInt(clique.get(j));
					if(v1!=v2) {
						edge += G.adjacency_matrix[v1][v2];		//increments when there is an edge between vertex v1 and v2
						numEdge += G.adjacency_matrix[v1][v2];
						weight += G.edgevalues[v1][v2];
					}
				}
				if(min[0] > numEdge && i != clique.indexOf(node)) {
					min[0] = numEdge;
					min[1] = weight;
					min[2] = i;
				}
			}
			if(edge != n*(n-1)){
				clique.remove((int) min[2]);
				n--;
			} else {
				break;
			}
		}
		if(clique.size() < 3 || ! clique.contains(node))
			return ("null");
		else 
			return clique.toString();
	}
	void removeRedundant(Hashtable h) {
		for(Enumeration<LinkedList> e=h.elements(); e.hasMoreElements();) {
			LinkedList l = e.nextElement();
			while(l.contains("null")) {
				l.remove("null");
			}
			for(int i=0; i< l.size()-1; i++) {
				int j=i+1;
				while(j < l.size()) {
					if(l.get(i).equals(l.get(j)))
						l.remove(j);
					else
						j++;
				}
			}
		}
	}
	
	String mergecliques(String c1, String c2) {
		c1 = c1.substring(1,c1.length()-1);
		c2 = c2.substring(1,c2.length()-1);
		c1 = c1 + ", " + c2;							//combines list
		//List<String> c = new LinkedList<String>(Arrays.asList(c1.split(", ")));
		List<String> c = String2List(c1);
		for(int i=0; i< c.size()-1; i++) {
			int j=i+1;
			while(j<c.size()) {
				if(c.get(i).equals(c.get(j)))
					c.remove(j);				//removes duplicate nodes 
				else 
					j++;
			}
		}
		Object[] temp = c.toArray();
		Arrays.sort(temp);						//sorts the nodes in order
		String str = Arrays.toString(temp);
		return str;
	}
	
	boolean quasiclique(String clique1, String clique2) {
		double rab, RAB;
		int intersect=0, min, max;
		//List<String> c1 = new LinkedList<String>(Arrays.asList(clique1.substring(1,clique1.length()-1).split(", ")));
		//List<String> c2 = new LinkedList<String>(Arrays.asList(clique2.substring(1,clique2.length()-1).split(", ")));
		List<String> c1 = String2List(clique1);
		List<String> c2 = String2List(clique2);
		
		if(c1.size() >= c2.size()) {
			max = c1.size();
			min = c2.size();
		} else {
			max = c2.size();
			min = c1.size();
		}
		for(int i=0; i<c1.size(); i++) {
			for(int j=0; j<c2.size(); j++) {
				if(c1.get(i).equals(c2.get(j))) {
					intersect++;
					continue;
				}
			}
		}
		rab = (double) intersect/min;
		RAB = (double) intersect/max;
		if(rab >.9 && RAB > .7) 
			return true;
		else 
			return false;
		
	}
	
	double clusterScore(weighed_matrix best, LinkedList<weighed_matrix> cluster) {
		int n = best.motif.length;	//number of sequences in the best motif
		int N = 0;					//number of sequences in cluster
		int L = best.L;				//L-length motif in best motif
		double x=0;
		double clusterscore = 0;
		for(int i=0; i<cluster.size(); i++) {
			N += cluster.get(i).motif.length;
		}
		
		for(int i=0; i<L; i++) 
			x += best.I(i)/L;
		
		clusterscore = (n-Math.log(N))*Math.exp(x);
		
		return clusterscore;
		
	}
	
	List<String> String2List(String s) {
		return (new LinkedList<String>(Arrays.asList(s.substring(1,s.length()-1).split(", "))));
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		minimcl m = new minimcl();
		SparseMatrix SM;
		predictMotif inputmotifs = new predictMotif();				//object
		File file = new File("motif.txt");							//input file
		Graph simGraph;												//graph of the similarity for pairwise comparisons of motifs
		ArrayList<String[]> listMotifs = new ArrayList<String[]>();	//array list of arrays where each node in the list is an array of the sequences found for that motif
		
		listMotifs = inputmotifs.fileIN(file, listMotifs);
		weighed_matrix[] matrixController = new weighed_matrix[listMotifs.size()];
		for(int i=0; i< listMotifs.size(); i++) 
			matrixController[i]= new weighed_matrix(listMotifs.get(i));					//motifs are put into a controller array 
		simGraph = new Graph(listMotifs.size());										//creates a graph of the correct size
		simGraph = inputmotifs.simgraph(matrixController);								//adds edges
		
		/*/artifical edges added for testing purpose
		simGraph.addedge(1,3);
		simGraph.addedge(3,1);
		simGraph.addedge(2,4);
		simGraph.addedge(4,2);
		*****************/
		
/*		m.addloops(simGraph.adjacency_matrix);											//start setting up for mcl by adding loops to graph nodes
		simGraph.adjacency_matrix = m.make_stochastic(simGraph.adjacency_matrix);
		SM = new SparseMatrix(simGraph.adjacency_matrix);*/
		m.addloops(simGraph.edgevalues);
		simGraph.edgevalues = m.make_stochastic(simGraph.edgevalues);
		SM = new SparseMatrix(simGraph.edgevalues);
		
		SM = m.mcl(SM, 2);
		String[] clusters = m.getClusters(SM);
		

		
		//Hashtable<String, String[]> cliquehash = new Hashtable<String, String[]>();
		Hashtable<String, LinkedList<String>> cliquehash = new Hashtable<String, LinkedList<String>>();
		for(int i=0; i<clusters.length; i++) {		//create a hashtable of the subgraph and cliques which are stored as strings of the nodes
			if( clusters[i] != ""){
				System.out.println(clusters[i]);
				List<String> c = new LinkedList<String>(Arrays.asList(clusters[i].split(" ")));
				int n = c.size();
				LinkedList<String> cliques = new LinkedList<String>();
				cliquehash.put(clusters[i], cliques);
				for(int j=0; j<n; j++) {
					 c = new LinkedList<String>(Arrays.asList(clusters[i].split(" ")));
					//cliquehash.get(clusters[i])[j] = inputmotifs.checkclique(c, simGraph);
					 cliquehash.get(clusters[i]).add(j, inputmotifs.checkclique(c.get(j), c, simGraph));	//cliques are added according to their subgraph. key = subgraph, value = cliques in subgraph
				}
				
				
			}
		}
		inputmotifs.removeRedundant(cliquehash);
		for(Enumeration<LinkedList<String>> e=cliquehash.elements(); e.hasMoreElements();) {
			LinkedList<String> l = e.nextElement();
			for(int i=0; i<l.size()-1; i++){
				int j=i+1;
				while(j<l.size()) {
					String clique1 = l.get(i);
					String clique2 = l.remove(j);
					if(inputmotifs.quasiclique(clique1, clique2))				//merges cliques after seeing if there is enough overlap
						l.set(i, inputmotifs.mergecliques(clique1, clique2));	//cliques are merged to the first clique in the list order
					else {
						l.add(j, clique2);										//if cliques are not merged, then cliques are left in original positions
						j++;
					}
				}
			}
		}
		
		//finds cluster score for a clique. future versions will find it for the cluster as a whole
		for(Enumeration<LinkedList<String>> e=cliquehash.elements(); e.hasMoreElements();) {
			LinkedList<String> cliques = e.nextElement();
			while(cliques.size() > 0){
				List<String> l = inputmotifs.String2List(cliques.pop());
				System.out.print("for clique " + l + " clusterscore: ");
				LinkedList<weighed_matrix> cluster = new LinkedList<weighed_matrix>();
				int best = Integer.parseInt(l.remove(0));
				while(l.size() > 0) {
					cluster.add(matrixController[Integer.parseInt(l.remove(0))]);
				}
				System.out.println(inputmotifs.clusterScore(matrixController[best], cluster));	//generate a clusterscore for each clique found 
			}
		}
		
		System.out.println("done");
	}

}
