int min_representation(char *s, int n) {
	for (int i=0; i<n; ++i) {
        s[i+n] = s[i];
    }
    s[n*2] = '\0';

	int i = 0, j = 1, k = 0;
    while(i < n && j < n) {
		k = 0;
		while(k < n && s[i+k] == s[j+k]) k++;
		if(k >= n) break;
		if(s[i+k] > s[j+k]) {
            i = std:: max (i + k + 1, j + 1);
        } else {
            j = std:: max (i + 1, j + k + 1);
        }
	}
	return std::min(i ,j);
}
