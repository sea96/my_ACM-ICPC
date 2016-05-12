public class ReadDemo {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		InputStream inputstream = System.in;
		OutputStream outputstream = System.out;
		InputReader in = new InputReader (inputstream);
		PrintWriter out = new PrintWriter (outputstream);
		//
		out.close ();
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
    }
}
