import java.io.*;
import java.util.*;
import java.math.*;

public class ReadDemo  {
	public static void main(String[] args) {
		//Scanner in = new Scanner(new BufferedInputStream(System.in));
		
		InputStream inputstream = System.in;
		OutputStream outputstream = System.out;
		
		InputReader in = new InputReader(inputstream);
		PrintWriter out = new PrintWriter(outputstream);
		
		Task solver = new Task();
		solver.solve(1, in, out);
		
		out.close ();
	}
	
	static class Task {
		public void solve(int testNumber, InputReader in, PrintWriter out) {
			// solve here
		}
	}
	
	static class InputReader {
        public BufferedReader reader;
        public StringTokenizer tokenizer;
        
        public InputReader(InputStream stream) {
            reader = new BufferedReader(new InputStreamReader(stream), 32768);
            tokenizer = null;
        }
        
        public String next() {
            while (tokenizer == null || !tokenizer.hasMoreTokens()) {
                try {
                    tokenizer = new StringTokenizer(reader.readLine());
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
            return tokenizer.nextToken();
        }
        
        public int nextInt() {
            return Integer.parseInt(next());
        }
        
        public BigInteger nextBigInteger() {
        	return new BigInteger(next());
        }
    }
}
