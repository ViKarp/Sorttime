#include <iostream>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
typedef std::chrono::high_resolution_clock Clock;
void hibbardsort(auto x[], auto n)       //SHELL SORT - HIBBARD'S INCREMENTS
{
    int i, j,k, increment, temp;
    int swp=0,comp=0;
    int val;
    val=(int)log(n+1)/log(2);
    increment =(int)pow(2,val)-1;
    // increment =(int)9* pow(4,n)- 9* pow(2,n);
    while (increment > 0)
    {
        for (i=0;i<increment;i++)
        {
            for(j=0;j<n;j+=increment)
            {
                temp=x[j];
                for(k=j-increment;k>=0&&temp<x[k];k-=increment)
                {
                    comp++;
                    swp++;
                    x[k+increment]=x[k];
                }
                x[k+increment]=temp;
                swp++;
            }
        }
        comp++;
        val--;
        if(increment!=1)
            increment=pow(2,val)-1;
        else
            increment = 0;
    }
}
void rasheska(auto *array, auto size )
{
    if (!array || !size)
        return;

    int jump = size;
    bool swapped = true;

    while (jump > 1 || swapped)
    {
        if (jump > 1)
            jump = (int)(jump / 1.25);
        swapped = false;
        for (int i = 0; i + jump < size; i++)
            if (array[i] > array[i + jump])
                std::swap(array[i],array[i + jump]), swapped = true;
    }
}
int getRN(int min, auto max)
{
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
    // Равномерно распределяем рандомное число в нашем диапазоне
    return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}
void merge(auto a[], auto fsize)
{
    if (fsize < 2)return;
    merge(a, fsize / 2);
    merge(&a[fsize / 2], fsize - (fsize / 2));
    int* buf = new int[fsize];
    int idbuf = 0, idl = 0, idr = fsize / 2 ;
    while ((idl < fsize / 2) && (idr < fsize))
        if (a[idl] < a[idr])
            buf[idbuf++] = a[idl++];
        else
            buf[idbuf++] = a[idr++];
    while (idl < fsize / 2) buf[idbuf++] = a[idl++];
    while (idr < fsize) buf[idbuf++] = a[idr++];
    for (idl = 0; idl < fsize; idl++)a[idl] = buf[idl];
    delete[]buf;
}
void quicksort(auto x[],auto first,auto last){ //Разбиение Хоара
    int pivot,j,temp,i;

    if(first<last){
        pivot=first;
        i=first;
        j=last;

        while(i<j){
            while(x[i]<=x[pivot]&&i<last)
                i++;
            while(x[j]>x[pivot])
                j--;
            if(i<j){
                temp=x[i];
                x[i]=x[j];
                x[j]=temp;
            }
        }

        temp=x[pivot];
        x[pivot]=x[j];
        x[j]=temp;
        quicksort(x,first,j-1);
        quicksort(x,j+1,last);

    }
}
void siftDown(auto *numbers, auto root, auto bottom)
{
    int maxChild; // индекс максимального потомка
    int done = 0; // флаг того, что куча сформирована
    // Пока не дошли до последнего ряда
    while ((root * 2 <= bottom) && (!done))
    {
        if (root * 2 == bottom)    // если мы в последнем ряду,
            maxChild = root * 2;    // запоминаем левый потомок
            // иначе запоминаем больший потомок из двух
        else if (numbers[root * 2] > numbers[root * 2 + 1])
            maxChild = root * 2;
        else
            maxChild = root * 2 + 1;
        // если элемент вершины меньше максимального потомка
        if (numbers[root] < numbers[maxChild])
        {
            int temp = numbers[root]; // меняем их местами
            numbers[root] = numbers[maxChild];
            numbers[maxChild] = temp;
            root = maxChild;
        }
        else // иначе
            done = 1; // пирамида сформирована
    }
}
void heapSort(auto *numbers, auto array_size)
{
    // Формируем нижний ряд пирамиды
    for (int i = (array_size / 2); i >= 0; i--)
        siftDown(numbers, i, array_size - 1);
    // Просеиваем через пирамиду остальные элементы
    for (int i = array_size - 1; i >= 1; i--)
    {
        int temp = numbers[0];
        numbers[0] = numbers[i];
        numbers[i] = temp;
        siftDown(numbers, 0, i - 1);
    }
}
void selectionSort(auto *num, auto size)
{
    int min, temp; // для поиска минимального элемента и для обмена
    for (int i = 0; i < size - 1; i++)
    {
        min = i; // запоминаем индекс текущего элемента
        // ищем минимальный элемент чтобы поместить на место i-ого
        for (int j = i + 1; j < size; j++)  // для остальных элементов после i-ого
        {
            if (num[j] < num[min]) // если элемент меньше минимального,
                min = j;       // запоминаем его индекс в min
        }
        temp = num[i];      // меняем местами i-ый и минимальный элементы
        num[i] = num[min];
        num[min] = temp;
    }
}
bool check(auto x[],int n) {
    bool Flag = true;
    for (int i = 1; i<n; i++) {
        if (x[i]<x[i-1]) {
            //std::cout << x[i-1] << " "<< x[i] << "\n";
            Flag = false;
        }
    }
    return Flag;
}
void insertionsort(int a[], int low, int high)
{
    // start from the second element in the subarray
    // (the element at index `low` is already sorted)
    for (int i = low + 1; i <= high; i++)
    {
        int value = a[i];
        int j = i;

        // find index `j` within the sorted subset a[0…i-1]
        // where element a[i] belongs
        while (j > low && a[j - 1] > value)
        {
            a[j] = a[j - 1];
            j--;
        }

        // Note that the subarray `a[j…i-1]` is shifted to
        // the right by one position, i.e., `a[j+1…i]`

        a[j] = value;
    }
}

// Function to partition the array using Lomuto partition scheme
int partition(int a[], int low, int high)
{
    // Pick the rightmost element as a pivot from the array
    int pivot = a[high];

    // elements less than the pivot will be pushed to the left of `pIndex`
    // elements more than the pivot will be pushed to the right of `pIndex`
    // equal elements can go either way
    int pIndex = low;

    // each time we find an element less than or equal to the pivot, `pIndex`
    // is incremented, and that element would be placed before the pivot.
    for (int i = low; i < high; i++)
    {
        if (a[i] <= pivot)
        {
            std::swap(a[i], a[pIndex]);
            pIndex++;
        }
    }

    // swap `pIndex` with pivot
    std::swap (a[pIndex], a[high]);

    // return `pIndex` (index of the pivot element)
    return pIndex;
}

// Quicksort randomized partition to rearrange elements across pivot
int randPartition(int a[], int low, int high)
{
    // choose a random index between `[low, high]`
    int pivotIndex = rand() % (high - low + 1) + low;

    // swap the end element with the element present at a random index
    std::swap(a[pivotIndex], a[high]);

    // call the partition procedure
    return partition(a, low, high);
}

// Function to perform heapsort on the given range of elements
void heapsort1(int *begin, int *end)
{
    std::make_heap(begin, end);
    std::sort_heap(begin, end);
}

// Function to perform introsort on the given array
void introsort(int a[], int *begin, int *end, int maxdepth)
{
    // perform insertion sort if partition size is 16 or smaller
    if ((end - begin) < 16) {
        insertionsort(a, begin - a, end - a);
    }
        // perform heapsort if the maximum depth is 0
    else if (maxdepth == 0) {
        heapsort1(begin, end + 1);
    }
    else {
        // otherwise, perform Quicksort
        int pivot = randPartition(a, begin - a, end - a);
        introsort(a, begin, a + pivot - 1, maxdepth - 1);
        introsort(a, a + pivot + 1, end, maxdepth - 1);
    }
}

int main() {
    std::ofstream outf("sorttime.txt");
    std::ofstream outf1("sorttime2.txt");
    std::ofstream outf2("sorttime3.txt");
    std::ofstream outf3("sorttime4.txt");
    std::ofstream outf4("sorttime5.txt");
    for (int n = 3000; n <= 100000; n+=1000) {
        std::string text = std::to_string(n);
        std::vector<std::vector<double>> b(7);
        for (int j = 0; j<10; j++) {
            int *arr = new int[n];                  //создание массива из n элементов
            int *brr = new int[n];
            int *crr = new int[n];
            int *drr = new int[n];
            int *err = new int[n];
            int *frr = new int[n];
            int *hrr = new int[n];
            for (int i = 0; i < n; i++) {
                arr[i] = rand();
                brr[i] = arr[i];
                crr[i] = arr[i];
                drr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }
            auto t1 = Clock::now();
            quicksort(arr, 0, n-1);
            auto t2 = Clock::now();
            b[0].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            merge(frr, n);
            t2 = Clock::now();
            b[1].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            rasheska(brr, n);
            t2 = Clock::now();
            b[2].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            hibbardsort(crr, n);
            t2 = Clock::now();
            b[3].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            selectionSort(drr, n);
            t2 = Clock::now();
            b[4].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            heapSort(err, n);
            t2 = Clock::now();
            b[5].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            introsort(hrr,hrr,hrr+n-1,log(n) * 2);
            t2 = Clock::now();
            b[6].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
        }
        for (int i = 0; i<7;i++) {
            std::vector<double> k = b[i];
            std::sort(k.begin(), k.end());
            text+=" "+std::to_string(round( (k[5] + k[4] ) / 2.0*1000)/1000);
        }
        text+="\n";
        outf << text << std::endl;
    }
    for (int n = 200000; n <= 5000000; n+=100000) {
        std::string text = std::to_string(n);
        std::vector<std::vector<double>> b(7);
        for (int j = 0; j<10; j++) {
            int *arr = new int[n];                  //создание массива из n элементов
            int *brr = new int[n];
            int *err = new int[n];
            int *frr = new int[n];
            int *hrr = new int[n];
            for (int i = 0; i < n; i++) {
                arr[i] = rand();
                brr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }

            auto t1 = Clock::now();
            quicksort(arr, 0, n-1);
            auto t2 = Clock::now();
            b[0].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            merge(frr, n);
            t2 = Clock::now();
            b[1].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            rasheska(brr, n);
            t2 = Clock::now();
            b[2].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            b[3].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            b[4].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            heapSort(err, n);
            t2 = Clock::now();
            b[5].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            introsort(hrr,hrr,hrr+n-1,log(n) * 2);
            t2 = Clock::now();
            b[6].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
        }
        for (int i = 0; i<7;i++) {
            std::vector<double> k = b[i];
            std::sort(k.begin(), k.end());
            text+=" "+std::to_string(round( (k[5] + k[4] ) / 2.0*1000)/1000);
        }
        text+="\n";
        outf << text << std::endl;
    }
    for (int n = 10000000; n <= 20000000; n+=1000000) {
        std::string text = std::to_string(n);
        std::vector<std::vector<double>> b(7);
        for (int j = 0; j<10; j++) {
            int *arr = new int[n];                  //создание массива из n элементов
            int *brr = new int[n];
            int *err = new int[n];
            int *frr = new int[n];
            int *hrr = new int[n];
            for (int i = 0; i < n; i++) {
                arr[i] = rand();
                brr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }
            auto t1 = Clock::now();
            quicksort(arr, 0, n-1);
            auto t2 = Clock::now();
            b[0].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            merge(frr, n);
            t2 = Clock::now();
            b[1].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            rasheska(brr, n);
            t2 = Clock::now();
            b[2].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            heapSort(err, n);
            t2 = Clock::now();
            b[3].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            b[4].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            b[5].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            introsort(hrr,hrr,hrr+n-1,log(n) * 2);
            t2 = Clock::now();
            b[6].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
        }
        for (int i = 0; i<7;i++) {
            std::vector<double> k = b[i];
            std::sort(k.begin(), k.end());
            text+=" "+std::to_string(round( (k[5] + k[4] ) / 2.0*1000)/1000);
        }

        text+="\n";
        outf << text << std::endl;
    }
    //Второй круг. С почти отсортированном массивом
    for (int n = 30000; n <= 100000; n+=1000) {
        std::string text = std::to_string(n);
        std::vector<std::vector<double>> b(7);
        for (int j = 0; j<10; j++) {
            int *arr = new int[n];                  //создание массива из n элементов
            int *brr = new int[n];
            int *crr = new int[n];
            int *drr = new int[n];
            int *err = new int[n];
            int *frr = new int[n];
            int *hrr = new int[n];
            arr[0] = RAND_MAX;
            brr[0] = arr[0];
            crr[0] = arr[0];
            drr[0] = arr[0];
            err[0] = arr[0];
            frr[0] = arr[0];
            hrr[0] = arr[0];
            for (int i = 1; i <n/2; i++) {
                arr[i] = RAND_MAX;
                brr[i] = arr[i];
                crr[i] = arr[i];
                drr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }

            for (int i = n/2; i <n-1000; i++) {
                arr[i] = rand();
                brr[i] = arr[i];
                crr[i] = arr[i];
                drr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }
            for (int i = n-1000; i <n; i++) {
                arr[i] = RAND_MAX;
                brr[i] = arr[i];
                crr[i] = arr[i];
                drr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }
            auto t1 = Clock::now();
            quicksort(arr, 0, n-1);
            auto t2 = Clock::now();
            b[0].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            merge(frr, n);
            t2 = Clock::now();
            b[1].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            rasheska(brr, n);
            t2 = Clock::now();
            b[2].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            hibbardsort(crr, n);
            t2 = Clock::now();
            b[3].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            selectionSort(drr, n);
            t2 = Clock::now();
            b[4].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            heapSort(err, n);
            t2 = Clock::now();
            b[5].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            introsort(hrr,hrr,hrr+n-1,log(n) * 2);
            t2 = Clock::now();
            b[6].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
        }
        for (int i = 0; i<7;i++) {
            std::vector<double> k = b[i];
            std::sort(k.begin(), k.end());
            text+=" "+std::to_string(round( (k[5] + k[4] ) / 2.0*1000)/1000);
        }
        text+="\n";
        outf1 << text << std::endl;
    }
    for (int n = 30000; n <= 100000; n+=1000) {
        std::string text = std::to_string(n);
        std::vector<std::vector<double>> b(7);
        for (int j = 0; j<10; j++) {
            int *arr = new int[n];                  //создание массива из n элементов
            int *brr = new int[n];
            int *crr = new int[n];
            int *drr = new int[n];
            int *err = new int[n];
            int *frr = new int[n];
            int *hrr = new int[n];
            arr[0] = RAND_MAX;
            brr[0] = arr[0];
            crr[0] = arr[0];
            drr[0] = arr[0];
            err[0] = arr[0];
            frr[0] = arr[0];
            hrr[0] = arr[0];
            for (int i = 1; i < n; i++) {
                arr[i] = getRN(arr[i-1]-1000,arr[i-1]);
                brr[i] = arr[i];
                crr[i] = arr[i];
                drr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }
            auto t1 = Clock::now();
            quicksort(arr, 0, n-1);
            auto t2 = Clock::now();
            b[0].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            merge(frr, n);
            t2 = Clock::now();
            b[1].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            rasheska(brr, n);
            t2 = Clock::now();
            b[2].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            hibbardsort(crr, n);
            t2 = Clock::now();
            b[3].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            selectionSort(drr, n);
            t2 = Clock::now();
            b[4].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            heapSort(err, n);
            t2 = Clock::now();
            b[5].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            introsort(hrr,hrr,hrr+n-1,log(n) * 2);
            t2 = Clock::now();
            b[6].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
        }
        for (int i = 0; i<7;i++) {
            std::vector<double> k = b[i];
            std::sort(k.begin(), k.end());
            text+=" "+std::to_string(round( (k[5] + k[4] ) / 2.0*1000)/1000);
        }
        text+="\n";
        outf3 << text << std::endl;
    }
    for (int n = 30000; n <= 100000; n+=1000) {
        std::string text = std::to_string(n);
        std::vector<std::vector<double>> b(7);
        for (int j = 0; j<10; j++) {
            int *arr = new int[n];                  //создание массива из n элементов
            int *brr = new int[n];
            int *crr = new int[n];
            int *drr = new int[n];
            int *err = new int[n];
            int *frr = new int[n];
            int *hrr = new int[n];
            arr[0] = RAND_MAX;
            brr[0] = arr[0];
            crr[0] = arr[0];
            drr[0] = arr[0];
            err[0] = arr[0];
            frr[0] = arr[0];
            hrr[0] = arr[0];
            for (int i = 0; i < n/1.5; i++) {
                arr[i] = getRN(arr[i-1]-1000,arr[i-1]);
                brr[i] = arr[i];
                crr[i] = arr[i];
                drr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }
            for (int i = n/1.5; i < n; i++) {
                arr[i] = rand();
                brr[i] = arr[i];
                crr[i] = arr[i];
                drr[i] = arr[i];
                err[i] = arr[i];
                frr[i] = arr[i];
                hrr[i] = arr[i];
            }
            auto t1 = Clock::now();
            quicksort(arr, 0, n-1);
            auto t2 = Clock::now();
            b[0].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            merge(frr, n);
            t2 = Clock::now();
            b[1].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            rasheska(brr, n);
            t2 = Clock::now();
            b[2].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            hibbardsort(crr, n);
            t2 = Clock::now();
            b[3].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            selectionSort(drr, n);
            t2 = Clock::now();
            b[4].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            heapSort(err, n);
            t2 = Clock::now();
            b[5].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
            t1 = Clock::now();
            introsort(hrr,hrr,hrr+n-1,log(n) * 2);
            t2 = Clock::now();
            b[6].push_back(std::chrono::round<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
        }
        for (int i = 0; i<7;i++) {
            std::vector<double> k = b[i];
            std::sort(k.begin(), k.end());
            text+=" "+std::to_string(round( (k[5] + k[4] ) / 2.0*1000)/1000);
        }
        text+="\n";
        outf4 << text << std::endl;
    }

    outf.close();
    outf1.close();
    return 0;
}