# Fast Fourier Transform (FFT) Cooley-Tukey Radix-2

Tugas Analisis Algoritma Komputasi (Magister Teknik Informatika, Universitas Hasanuddin, Semester Ganjil 2026)
Algoritma No. 4: **Fast Fourier Transform**

| | |
|---|---|
| Nama | Besse Zahra Nurani A |
| NIM | D082261030 |
| Program Studi | Magister Teknik Informatika, Universitas Hasanuddin |
| Laporan | *Analisis Algoritma Fast Fourier Transform Cooley-Tukey Radix-2: Penurunan Kompleksitas, Pembuktian Kebenaran, dan Evaluasi Empiris* |

## Deskripsi

Repositori ini berisi implementasi dari nol (tanpa memanggil fungsi FFT pustaka) algoritma FFT Cooley-Tukey radix-2
*decimation-in-time* versi rekursif dengan paradigma *divide and conquer*. FFT menghitung Transformasi Fourier Diskret (DFT)

```
X_k = sum_{j=0}^{n-1} x_j * w_n^(jk),   w_n = e^(-2*pi*i/n),   k = 0, 1, ..., n-1
```

dalam waktu **Θ(n log n)**, dibandingkan dengan **Θ(n²)** pada DFT naif. Analisis teori (contoh manual, pseudocode,
bukti kebenaran dengan induksi, rekurensi `T(n) = 2T(n/2) + Θ(n)`, kompleksitas ruang Θ(n)) ada di laporan;
repositori ini memuat kode dan hasil eksperimennya.

## Isi Repositori

| Berkas | Keterangan |
|---|---|
| `fft.py` | Program utama: `dft_naive()`, `fft()`, `ifft()`, enam kasus uji, studi kasus, eksperimen waktu, dan pembuatan grafik |
| `eksperimen_tambahan.py` | Eksperimen opsional: kasus tepi (n = 1, 2, 3) dan pengaruh jenis serta tipe data terhadap waktu FFT |
| `hasil_fft.png` | Grafik hasil eksperimen waktu (skala log-log) terhadap kurva teori n log₂ n dan n² |
| `README.md` | Dokumen ini |

## Prasyarat

Pustaka yang dibutuhkan:

| Pustaka | Versi | Kegunaan |
|---|---|---|
| Python | 3.12 (diuji pada 3.12.3) | Bahasa pemrograman |
| numpy | 2.4.4 (diuji); versi 1.26 atau lebih baru seharusnya cukup | Larik, statistik, dan pembanding `numpy.fft` |
| matplotlib | 3.8 atau lebih baru | Pembuatan grafik `hasil_fft.png` |

Modul `cmath`, `math`, dan `time` adalah pustaka standar Python sehingga tidak perlu dipasang.

## Cara Menjalankan

```bash
# 1. (opsional) buat virtual environment
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 2. pasang pustaka
pip install numpy matplotlib

# 3. jalankan program utama
python fft.py

# 4. (opsional) jalankan eksperimen tambahan
python eksperimen_tambahan.py
```

Catatan: `eksperimen_tambahan.py` mengimpor `fft()` dari `fft.py`, jadi kedua berkas harus berada di folder yang sama.
Eksperimen waktu pada `fft.py` memakan waktu beberapa menit karena DFT naif diukur sampai n = 2048 (sekitar 5 detik per jalan).

## Keluaran Program

`python fft.py` mencetak dan menyimpan:

1. **Kasus uji (6 kasus, semuanya lulus):**
   - x = [1, 2, 3, 4] menghasilkan (10, -2+2i, -2, -2-2i)
   - impuls n = 8 menghasilkan semua keluaran bernilai 1
   - konstan [1, 1, 1, 1] menghasilkan (4, 0, 0, 0)
   - acak n = 16 dibandingkan dengan DFT naif (selisih < 10⁻⁹)
   - acak n = 16 dibandingkan dengan `numpy.fft`
   - `ifft(fft(x))` kembali ke x semula
2. **Studi kasus:**
   - Perkalian polinomial dengan FFT: (1 + 2x + 3x²)(4 + 5x) = 4 + 13x + 22x² + 15x³ (sama dengan perkalian biasa)
   - Deteksi frekuensi sinyal sin(2π·50t) + 0,5·sin(2π·120t), disampel 1.024 Hz sebanyak 1.024 titik, menghasilkan dua puncak pada 50 Hz dan 120 Hz
3. **Eksperimen waktu:** tabel waktu rata-rata ± simpangan baku dan rasio waktu/(n log₂ n)
4. **Grafik** `hasil_fft.png`

## Ringkasan Hasil Eksperimen

Ukuran input n = 2⁶ sampai 2¹⁵ (10 ukuran), 5 pengulangan per ukuran, waktu diukur dengan `time.perf_counter`.
Seed tetap: 42 (uji kebenaran) dan 2026 (eksperimen waktu). DFT naif hanya diukur sampai n = 2048.

| n | FFT sendiri (detik) | DFT naif (detik) | numpy.fft (detik) |
|---:|---:|---:|---:|
| 64 | 0,000304 | 0,005291 | 1,98e-05 |
| 128 | 0,000353 | 0,021176 | 1,56e-05 |
| 256 | 0,000814 | 0,106344 | 2,12e-05 |
| 512 | 0,002340 | 0,318646 | 2,08e-05 |
| 1024 | 0,003209 | 1,229890 | 3,75e-05 |
| 2048 | 0,007205 | 5,406277 | 4,82e-05 |
| 4096 | 0,015159 | - | 6,15e-04 |
| 8192 | 0,031854 | - | 5,49e-04 |
| 16384 | 0,066489 | - | 1,16e-03 |
| 32768 | 0,145329 | - | 2,90e-03 |

Temuan utama:

- Kemiringan regresi log-log sekitar **1,02** untuk FFT dan **1,98** untuk DFT naif, sesuai teori Θ(n log n) dan Θ(n²).
- Rasio t/(n log₂ n) FFT relatif konstan (sekitar 3,0×10⁻⁷ untuk n ≥ 4096), yaitu konstanta tersembunyi pada Θ(n log n).
- Pada n = 2048, FFT buatan sendiri sekitar 750 kali lebih cepat daripada DFT naif.
- Pada n = 32768, implementasi Python murni sekitar 50 kali lebih lambat daripada `numpy.fft` (kode C teroptimasi).
- Galat absolut maksimum terhadap `numpy.fft` kecil: 1,4×10⁻¹⁵ (n = 64) sampai 1,0×10⁻¹³ (n = 16384).
- Waktu FFT tidak bergantung pada isi data (nol, konstan, impuls, acak, terurut), tetapi elemen bertipe `numpy.float64`
  membuat FFT sekitar 1,8 kali lebih lambat daripada `float` Python biasa.

Waktu eksekusi bergantung pada mesin, sehingga angka di mesin lain akan berbeda; yang dibandingkan adalah pola pertumbuhannya.

## Batasan

- Hanya FFT satu dimensi radix-2 *decimation-in-time* versi rekursif.
- Panjang input harus pangkat dua; selain itu program melempar `ValueError("Panjang input harus pangkat 2")`.
- Aritmetika *floating-point* presisi ganda; memori Θ(n) pada versi rekursif.
- Pengukuran dilakukan pada satu mesin dan ukuran input maksimum 2¹⁵ karena keterbatasan kecepatan Python murni.

## Referensi Utama

- J. W. Cooley dan J. W. Tukey, "An algorithm for the machine calculation of complex Fourier series," *Mathematics of Computation*, 1965.
- T. H. Cormen dkk., *Introduction to Algorithms*, Bab 30 (Polynomials and the FFT).

Daftar referensi lengkap ada pada laporan.

## Pelaporan Penggunaan Alat Bantu AI

Penggunaan asisten AI dalam pengerjaan tugas ini dijelaskan pada Lampiran A laporan.
