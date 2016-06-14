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
