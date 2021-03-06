# parallelized-mst
Parallelize the Minimum Spanning Tree Algorithm with using OpenMP and OpenMPI

Kelompok 35:<br>
- 13518044 Jun Ho Choi Hedyatmo
- 13518149 Stefanus Gusega Gunawan

Cara Kerja Algoritma OpenMP Untuk MST:<br>
disini algoritma yang digunakan adalah algoritma prim, yang versi O(v^2). Untuk referensinya dari cp-algorithms.com. Intinya adalah dalam v step, kita harus mengambil minimum weight yang terdapat pada vertex yang belum divisit, jadi esensinya adalah mencari global minimum dari suatu array dan indeksnya sebanyak v kali. Untuk mencari global minimum dari suatu array dapat diparalelisasikan karena minimum dari suatu subarray dapat dicari secara independen, dan dapat dibandingkan langsung dengan minimum subarray lainnya untuk mencari global minimum, karena hal inilah, paralelisasi mungkin dilakukan. 
<br>
Cara Kerja Algoritma OpenMPI untuk MST:<br>
Algoritma yang digunakan adalah Kruskal, dengan referensi dari berbagai sumber seperti github.com, mpitutorial.com, dan geeksforgeeks.org. Intinya kegunaan OpenMPI di sini adalah berguna untuk membagi list of edges dari graph dengan metode MPI_Scatter, lalu ketika dilakukan parallel sorting dengan menggunakan merge sort, antarproses akan saling berkirim-kiriman list of edges, hingga nantinya ketika sudah di-merge sort akan dimerge lagi. Alasan menggunakan merge sort adalah kompleksitasnya yang cukup ringan yaitu O(nlogn). Inti dari paralelisasi yang digunakan adalah, setiap edge list yang sudah terpartisi diurutkan di dalam proses-proses, lalu dimerge antarproses, dan hasilpun dikeluarkan.