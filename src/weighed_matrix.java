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
	String[] motif;
	Hashtable<String, double[]> basefreq;
	Hashtable<String, double[]> Pm;

	weighed_matrix() {
		basefreq 	= new Hashtable<String, double[]>();
		Pm 			= new Hashtable<String, double[]>();		
		initialize(basefreq);
		initialize(Pm);
	}
	
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
		hash.put("A", a);
		hash.put("C", c);
		hash.put("G", g);
		hash.put("T", t);
		return hash;
	}

	void freq() {
		for(int i=0; i< this.motif.length; i++) {
			for(int j=0; j<this.motif[i].length(); j++) {
				if(Pattern.matches("(A|a)", this.motif[i].substring(j,j+1))) {
					this.basefreq.get("A")[j] += 1;
				}else if (Pattern.matches("(C|c)",this.motif[i].substring(j,j+1))) {
					this.basefreq.get("C")[j] += 1;
				}else if (Pattern.matches("(G|g)", this.motif[i].substring(j,j+1))) {
					this.basefreq.get("G")[j] += 1;
				}else if (Pattern.matches("(T|t)", this.motif[i].substring(j,j+1))) {
					this.basefreq.get("T")[j] += 1;
				}
				this.basefreq.get("A")[j] += this.pseudocount/this.motif.length;
				this.basefreq.get("C")[j] += this.pseudocount/this.motif.length;
				this.basefreq.get("G")[j] += this.pseudocount/this.motif.length;
				this.basefreq.get("T")[j] += this.pseudocount/this.motif.length;
			}
		}
	}

	void buildPm() {
		String[] base = {"A","C","G","T"};
		for (int i =0; i<4; i++) {							//loop for each base
			int len = this.basefreq.get(base[i]).length;	//loops the for length L 
			for (int j=0; j<len;j++) {
				this.Pm.get(base[i])[j] = Math.log(Pm(base[i],j)/q(base[i]));	//Profile matrix
			}
		}
	}
	
	double Pm (String base, int i) {
		double a 	= this.basefreq.get("A")[i];
		double c 	= this.basefreq.get("C")[i];
		double g 	= this.basefreq.get("G")[i];
		double t 	= this.basefreq.get("T")[i];
		return (this.basefreq.get(base)[i]/(a+c+g+t)); //probability of base b at position i
	}
	
	double q (String base) {
		int len = this.basefreq.get(base).length;
		int a=0, c=0, g=0, t=0, nbase=0;
		for (int i=0; i<len; i++){
			a += this.basefreq.get("A")[i];
			c += this.basefreq.get("C")[i];
			g += this.basefreq.get("G")[i];
			t += this.basefreq.get("T")[i];
			nbase += this.basefreq.get(base)[i];
		}
		return (nbase/(a+c+g+t));
	}

	double I(int i){
		double I=0;
		Enumeration<String> b = Pm.keys();
		while( b.hasMoreElements()) {
			String base = b.nextElement(); 
			I = I + Pm.get(base)[i] * basefreq.get(base)[i];
		}
		return I;
	}
	
	double score(P, F) {
		
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
