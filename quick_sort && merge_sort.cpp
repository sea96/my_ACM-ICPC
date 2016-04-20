//quick_sort：快速排序，随机标杆，O(nlogn)
void quick_sort(int *a, int l, int r)   {
    if (l < r)  {
        swap (a[l+rand()%(r-l+1)], a[l]);
        int i = l, j = r, x = a[l];
        while (i < j)   {
            while (i < j && a[j] >= x)  j--;
            if (i < j)  a[i++] = a[j];
            while (i < j && a[i] < x)   i++;
            if (i < j)  a[j--] = a[i];
        }
        a[i] = x;
        quick_sort (a, l, i-1);
        quick_sort (a, i+1, r);
    }
}
//merge_sort：归并排序，O(nlogn)
void merge(int *a, int p, int q, int r)   {
    int n1 = q - p + 1, n2 = r - q;
    for (int i=1; i<=n1; ++i)   L[i] = a[p+i-1];
    for (int i=1; i<=n2; ++i)   R[i] = a[q+i];
    L[n1+1] = R[n2+1] = INF;
    for (int i=1, j=1, k=p; k<=r; ++k)  {   //求逆序数：cnt += n1 - i + 1;
        a[k] = (L[i] <= R[j]) ? L[i++] : R[j++];
    }
}
void merge_sort(int *a, int p, int r)   {
    if (p < r)  {
        int q = (p + r) >> 1;
        merge_Sort (a, p, q);
        merge_Sort (a, q+1, r);
        merge (a, p, q, r);
    }
}
/*
   srand (time (NULL));
   clock_t begin = clock ();
   quick_sort (a, 0, 10000000);
   clock_t end = clock ();
   printf ("%.2fMS\n", (double) (end - begin) * 1000 / CLOCKS_PER_SEC);
*/
