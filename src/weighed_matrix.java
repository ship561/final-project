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
	int L;
	String[] motif;
	Hashtable<String, double[]> basefreq;
	Hashtable<String, double[]> Prf;
	
	weighed_matrix(String[] motifs) {
		motif		= motifs;
		basefreq 	= new Hashtable<String, double[]>();		//base frequency
		Prf 		= new Hashtable<String, double[]>();		//profile matrix ln(prob(b,i)/q(b))
		L			= motifs[0].length();
		initialize(basefreq);
		initialize(Prf);
		freq();
		Prf();
	}

	Hashtable<String, double[]> initialize(Hashtable<String,double[]> hash) {
		double[] a = new double[L];
		double[] c = new double[L];
		double[] g = new double[L];
		double[] t = new double[L];
		hash.put("A", a);
		hash.put("C", c);
		hash.put("G", g);
		hash.put("T", t);
		return hash;
	}

	void freq() {	//frequency table
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
				this.basefreq.get("A")[j] += (double) this.pseudocount/this.motif.length;
				this.basefreq.get("C")[j] += (double) this.pseudocount/this.motif.length;
				this.basefreq.get("G")[j] += (double) this.pseudocount/this.motif.length;
				this.basefreq.get("T")[j] += (double) this.pseudocount/this.motif.length;
			}
		}
	}

	void Prf() {	//ln(pm(b,i)/q(b))
		String[] base = {"A","C","G","T"};
		for (int i =0; i<4; i++) {							//loop for each base
			int len = this.basefreq.get(base[i]).length;	//loops the for length L 
			for (int j=0; j<len;j++) {
				this.Prf.get(base[i])[j] = Math.log(pm(base[i],j)/q(base[i]));	//Profile matrix
			}
		}
	}
	
	double pm (String base, int i) {	//probability of base b at position i
		double p = this.basefreq.get(base)[i]/this.basefreq.get(base).length;
		return p;

	}
	
	double q (String base) {	//background frequency of base b
		int len = this.basefreq.get(base).length;
		int a=0, c=0, g=0, t=0, nbase=0;
		for (int i=0; i<len; i++){
			a += this.basefreq.get("A")[i];
			c += this.basefreq.get("C")[i];
			g += this.basefreq.get("G")[i];
			t += this.basefreq.get("T")[i];
			nbase += this.basefreq.get(base)[i];
		}
		return ((double) nbase/(a+c+g+t));
	}

	double I(int i){	//information stored in each column
		double I=0;
		Enumeration<String> b = Prf.keys();
		while( b.hasMoreElements()) {
			String base = b.nextElement(); 
			I = I + pm(base,i) * Prf.get(base)[i];
		}
		return I;
	}
	
}
