public void solve(int testNumber, InputReader in, PrintWriter out) {
	//求x=sqrt(n)，小数点后保留k位
	BigDecimal n = in.nextBigDecimal();
	int k = in.nextInt();
	
	BigDecimal x = BigDecimal.valueOf(Math.sqrt(n.doubleValue()));
	
	//Newton's method: x(n+1)=x(n)-f(x(n))/f'(x(n))
	for (int i=0; i<10; ++i) {
	x = x.add(n.divide(x, k+5, RoundingMode.HALF_UP));
	x = x.divide(BigDecimal.valueOf(2), k+5, RoundingMode.HALF_UP);
	}
	x = x.setScale(k, BigDecimal.ROUND_DOWN);
	out.println(x);
}
