import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.regex.Pattern;


public class weighed_matrix {
	int pseudocount = 2;

	void fileIN(File f) {           //reads file in
		//StringBuffer contents = new StringBuffer();
		BufferedReader reader = null;
		String[] input = new String[350];

		try {
			reader = new BufferedReader(new FileReader(f));
			String text = null;
			int j=0;
			while (reader.readLine() != null) {
				text = reader.readLine();
				if(! Pattern.matches("^>.*", text)) {
					input[j] = text;
					j++;
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		this.sequences=input;
	}
	Hashtable<String, double[]> initialize(Hashtable<String,double[]> hash) {
		double[] a = new double[wordlength];
		double[] c = new double[wordlength];
		double[] g = new double[wordlength];
		double[] t = new double[wordlength];
		double[] tot = new double[wordlength];
		hash.put("A", a);
		hash.put("C", c);
		hash.put("G", g);
		hash.put("T", t);
		hash.put("total", tot);
		return hash;
	}

	Hashtable<String, double[]> freq() {
		Hashtable<String,double[]>temp_table = new Hashtable<String,double[]>();
		initialize(temp_table); //keys = [A, C, G, T] array = pos[1-6]
		for(int i=0; i< this.sample.length; i++) {
			for(int j=0; j<this.sample[i].length(); j++) {
				if(Pattern.matches("(A|a)", this.sample[i].substring(j,j+1))) {
					temp_table.get("A")[j] += 1;
				}else if (Pattern.matches("(C|c)",this.sample[i].substring(j,j+1))) {
					temp_table.get("C")[j] += 1;
				}else if (Pattern.matches("(G|g)", this.sample[i].substring(j,j+1))) {
					temp_table.get("G")[j] += 1;
				}else if (Pattern.matches("(T|t)", this.sample[i].substring(j,j+1))) {
					temp_table.get("T")[j] += 1;
				}
				temp_table.get("A")[j] += this.pseudocount/this.sample.length;
				temp_table.get("C")[j] += this.pseudocount/this.sample.length;
				temp_table.get("G")[j] += this.pseudocount/this.sample.length;
				temp_table.get("T")[j] += this.pseudocount/this.sample.length;
			}
		}
		return temp_table;
	}

	Hashtable Pm (Hashtable<String, double[]> freq) {
		int a=0;
		int c=0;
		int g=0;
		int t=0;
		Hashtable<String, double[]> Pm = new Hashtable();
		initialize(Pm);
		Pm.remove("total");
		for(int i = 0; i<this.wordlength; i++) {
			a+=freq.get("A")[i];
			c+=freq.get("C")[i];
			g+=freq.get("G")[i];
			t+=freq.get("T")[i];
			freq.get("total")[i] += freq.get("A")[i] + freq.get("C")[i] + freq.get("G")[i] + freq.get("T")[i];
		}
		for(int i=0; i<this.wordlength; i++) {
			Pm.get("A")[i] = Math.log((freq.get("A")[i]/freq.get("total")[i])/((a-freq.get("A")[i])/(a+c+g+t)));
			Pm.get("C")[i] = Math.log((freq.get("C")[i]/freq.get("total")[i])/((c-freq.get("C")[i])/(a+c+g+t)));
			Pm.get("G")[i] = Math.log((freq.get("G")[i]/freq.get("total")[i])/((g-freq.get("G")[i])/(a+c+g+t)));
			Pm.get("T")[i] = Math.log((freq.get("T")[i]/freq.get("total")[i])/((t-freq.get("T")[i])/(a+c+g+t)));	
		}
		return Pm;
	}

	double I(int i, Hashtable<String, double[]> Pm, Hashtable<String, double[]> F){
		double I=0;
		Enumeration<String> b = Pm.keys();
		while( b.hasMoreElements()) {
			String base = b.nextElement(); 
			I = I + Pm.get(base)[i] * F.get(base)[i];
		}
		return I;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
