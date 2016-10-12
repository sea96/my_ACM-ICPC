int main() {
   srand(time(NULL));
   clock_t begin = clock();
   solve();
   clock_t end = clock();
   printf ("%.2fMS\n", (double)(end-begin)*1000/CLOCKS_PER_SEC);
   return 0;
}
