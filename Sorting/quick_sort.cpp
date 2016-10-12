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
