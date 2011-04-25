import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Enumeration;
import java.util.regex.Pattern;


public class predictMotif {

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
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
