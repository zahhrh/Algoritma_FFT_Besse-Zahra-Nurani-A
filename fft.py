"""
Fast Fourier Transform (FFT) - Cooley-Tukey Radix-2 (Divide & Conquer)
Tugas Analisis Algoritma Komputasi - Besse Zahra Nurani A (D082261030)

Cara menjalankan:
    pip install numpy matplotlib
    python fft.py

Isi file:
    1. dft_naive()   -> DFT langsung, O(n^2)  (pembanding & pembuktian kebenaran)
    2. fft()         -> FFT rekursif, O(n log n)  (INTI, ditulis sendiri)
    3. ifft()        -> FFT balik (dipakai untuk perkalian polinomial)
    4. uji_kasus()   -> 3+ kasus uji kecil yang bisa dicek manual
    5. studi_kasus() -> perkalian polinomial & deteksi frekuensi sinyal
    6. eksperimen()  -> ukur waktu, 5x ulang, rata-rata & simpangan baku, buat grafik

Sumber konsep: Cooley & Tukey (1965); Cormen dkk., Introduction to Algorithms, Bab 30.
numpy HANYA dipakai untuk array, plotting, dan sebagai pembanding (np.fft).
"""
import cmath          # cmath.exp untuk bilangan kompleks
import math
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")  # supaya grafik bisa disimpan tanpa jendela
import matplotlib.pyplot as plt


# ---------------------------------------------------------------
# 1. DFT NAIF: X[k] = sum_{j=0}^{n-1} x[j] * e^{-2*pi*i*j*k/n}
#    Dua loop bersarang -> O(n^2)
# ---------------------------------------------------------------
def dft_naive(x):
    n = len(x)
    hasil = []
    for k in range(n):                      # untuk tiap frekuensi k
        jumlah = 0
        for j in range(n):                  # jumlahkan semua sampel
            jumlah += x[j] * cmath.exp(-2j * math.pi * j * k / n)
        hasil.append(jumlah)
    return hasil


# ---------------------------------------------------------------
# 2. FFT REKURSIF (INTI)  - panjang input harus pangkat 2
#    Divide  : pisah indeks genap dan ganjil
#    Conquer : FFT kedua setengah secara rekursif
#    Combine : X[k]       = E[k] + w^k * O[k]
#              X[k + n/2] = E[k] - w^k * O[k]     (w = e^{-2*pi*i/n})
#    Rekurensi: T(n) = 2T(n/2) + Theta(n)  ->  Theta(n log n)
# ---------------------------------------------------------------
def fft(x):
    n = len(x)
    if n == 1:                              # basis rekursi
        return [x[0]]
    if n % 2 != 0:
        raise ValueError("Panjang input harus pangkat 2")

    genap = fft(x[0::2])                    # DIVIDE: indeks 0,2,4,...
    ganjil = fft(x[1::2])                   # DIVIDE: indeks 1,3,5,...

    hasil = [0] * n
    for k in range(n // 2):                 # COMBINE
        w = cmath.exp(-2j * math.pi * k / n)    # faktor putar (twiddle factor)
        hasil[k] = genap[k] + w * ganjil[k]
        hasil[k + n // 2] = genap[k] - w * ganjil[k]
    return hasil


# ---------------------------------------------------------------
# 3. FFT BALIK: ifft(X) = conj( fft( conj(X) ) ) / n
# ---------------------------------------------------------------
def ifft(X):
    n = len(X)
    konj = [v.conjugate() for v in X]
    y = fft(konj)
    return [v.conjugate() / n for v in y]


# ---------------------------------------------------------------
# Fungsi bantu
# ---------------------------------------------------------------
def hampir_sama(a, b, tol=1e-9):
    return all(abs(p - q) < tol for p, q in zip(a, b))


# ---------------------------------------------------------------
# 4. KASUS UJI KECIL (bisa dicek manual dengan kertas)
# ---------------------------------------------------------------
def uji_kasus():
    print("=== KASUS UJI ===")

    # Uji 1: n=4, x=[1,2,3,4]  ->  X=[10, -2+2i, -2, -2-2i]  (hitung manual)
    x = [1, 2, 3, 4]
    harapan = [10, -2 + 2j, -2, -2 - 2j]
    ok = hampir_sama(fft(x), harapan)
    print("Uji 1  x=[1,2,3,4]            :", "LULUS" if ok else "GAGAL")

    # Uji 2: impuls [1,0,0,0,0,0,0,0] -> semua keluaran bernilai 1
    x = [1, 0, 0, 0, 0, 0, 0, 0]
    ok = hampir_sama(fft(x), [1] * 8)
    print("Uji 2  impuls n=8             :", "LULUS" if ok else "GAGAL")

    # Uji 3: sinyal konstan [1,1,1,1] -> [4,0,0,0]
    x = [1, 1, 1, 1]
    ok = hampir_sama(fft(x), [4, 0, 0, 0])
    print("Uji 3  konstan n=4            :", "LULUS" if ok else "GAGAL")

    # Uji 4: acak n=16, bandingkan dengan DFT naif dan numpy
    rng = np.random.default_rng(42)         # seed tetap
    x = list(rng.random(16))
    ok1 = hampir_sama(fft(x), dft_naive(x))
    ok2 = hampir_sama(fft(x), list(np.fft.fft(x)))
    print("Uji 4  acak n=16 vs DFT naif  :", "LULUS" if ok1 else "GAGAL")
    print("Uji 5  acak n=16 vs np.fft    :", "LULUS" if ok2 else "GAGAL")

    # Uji 6: ifft(fft(x)) harus kembali ke x
    ok = hampir_sama(ifft(fft(x)), x)
    print("Uji 6  ifft(fft(x)) == x      :", "LULUS" if ok else "GAGAL")
    print()


# ---------------------------------------------------------------
# 5. STUDI KASUS
#    (a) Perkalian polinomial lewat FFT: O(n log n) vs O(n^2)
#    (b) Deteksi frekuensi dominan sinyal
# ---------------------------------------------------------------
def kali_polinomial_fft(a, b):
    n = 1
    while n < len(a) + len(b) - 1:          # cari pangkat 2 >= panjang hasil
        n *= 2
    A = fft(list(a) + [0] * (n - len(a)))   # evaluasi di akar satuan
    B = fft(list(b) + [0] * (n - len(b)))
    C = [p * q for p, q in zip(A, B)]       # perkalian titik demi titik
    c = ifft(C)                             # interpolasi
    return [round(v.real) for v in c[: len(a) + len(b) - 1]]


def kali_polinomial_naif(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i in range(len(a)):
        for j in range(len(b)):
            c[i + j] += a[i] * b[j]
    return c


def studi_kasus():
    print("=== STUDI KASUS ===")
    # (a) (1 + 2x + 3x^2) * (4 + 5x)  = 4 + 13x + 22x^2 + 15x^3
    a, b = [1, 2, 3], [4, 5]
    print("(a) Perkalian polinomial FFT  :", kali_polinomial_fft(a, b))
    print("    Perkalian polinomial naif :", kali_polinomial_naif(a, b))

    # (b) sinyal = sinus 50 Hz + sinus 120 Hz, laju sampling 1024 Hz, 1024 sampel
    fs, n = 1024, 1024
    t = [i / fs for i in range(n)]
    sinyal = [math.sin(2 * math.pi * 50 * ti) + 0.5 * math.sin(2 * math.pi * 120 * ti)
              for ti in t]
    spektrum = fft(sinyal)
    besar = [abs(v) for v in spektrum[: n // 2]]       # setengah pertama saja
    urut = sorted(range(len(besar)), key=lambda k: -besar[k])[:2]
    print("(b) Dua frekuensi dominan (Hz):", sorted(k * fs / n for k in urut))
    print()
    return t, sinyal, besar, fs, n


# ---------------------------------------------------------------
# 6. EKSPERIMEN WAKTU (>=5 ukuran, 5x ulang, rata-rata & std)
# ---------------------------------------------------------------
def ukur(fungsi, x, ulang=5):
    waktu = []
    for _ in range(ulang):
        mulai = time.perf_counter()
        fungsi(x)
        waktu.append(time.perf_counter() - mulai)
    return float(np.mean(waktu)), float(np.std(waktu, ddof=1))


def eksperimen():
    print("=== EKSPERIMEN WAKTU (detik, rata-rata +- std, 5 pengulangan) ===")
    rng = np.random.default_rng(2026)       # seed tetap
    ukuran = [2 ** p for p in range(6, 16)]  # 64 ... 32768
    batas_naif = 2 ** 11                    # DFT naif terlalu lambat di atas ini

    baris = []
    print(f"{'n':>7} | {'FFT sendiri':>20} | {'DFT naif':>20} | {'numpy.fft':>20}")
    for n in ukuran:
        x = list(rng.random(n))
        m_f, s_f = ukur(fft, x)
        m_n, s_n = ukur(dft_naive, x) if n <= batas_naif else (None, None)
        xa = np.array(x)
        m_p, s_p = ukur(np.fft.fft, xa)
        baris.append((n, m_f, s_f, m_n, s_n, m_p, s_p))
        naif_txt = f"{m_n:.6f} +- {s_n:.6f}" if m_n is not None else "-"
        print(f"{n:>7} | {m_f:.6f} +- {s_f:.6f} | {naif_txt:>20} | {m_p:.2e} +- {s_p:.1e}")
    print()
    return baris


def buat_grafik(baris, t, sinyal, besar, fs, n_sinyal):
    n = np.array([b[0] for b in baris], dtype=float)
    m_f = np.array([b[1] for b in baris])
    s_f = np.array([b[2] for b in baris])
    nn = [b for b in baris if b[3] is not None]
    n_n = np.array([b[0] for b in nn], dtype=float)
    m_n = np.array([b[3] for b in nn])

    # skala kurva teori agar berimpit di titik terbesar masing-masing
    teori_nlogn = n * np.log2(n)
    c1 = m_f[-1] / teori_nlogn[-1]
    teori_n2 = n_n ** 2
    c2 = m_n[-1] / teori_n2[-1]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].errorbar(n, m_f, yerr=s_f, fmt="o-", label="FFT (eksperimen)")
    ax[0].plot(n, c1 * teori_nlogn, "--", label=r"$c\,n\log_2 n$ (teori)")
    ax[0].plot(n_n, m_n, "s-", label="DFT naif (eksperimen)")
    ax[0].plot(n_n, c2 * teori_n2, ":", label=r"$c\,n^2$ (teori)")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("Ukuran input n"); ax[0].set_ylabel("Waktu (detik)")
    ax[0].legend(fontsize=8); ax[0].grid(True, which="both", alpha=0.3)
    ax[0].set_title("Waktu vs ukuran input")

    ax[1].plot(np.arange(len(besar)) * fs / n_sinyal, besar)
    ax[1].set_xlim(0, 200)
    ax[1].set_xlabel("Frekuensi (Hz)"); ax[1].set_ylabel("|X[k]|")
    ax[1].set_title("Spektrum: puncak di 50 Hz dan 120 Hz")
    ax[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("hasil_fft.png", dpi=200)
    print("Grafik disimpan: hasil_fft.png")

    # tampilkan rasio waktu/(n log n) untuk bukti pertumbuhan
    print("Rasio waktu/(n log2 n) FFT (harus kira-kira konstan):")
    for ni, mi in zip(n, m_f):
        print(f"  n={int(ni):>6}: {mi / (ni * math.log2(ni)):.3e}")


if __name__ == "__main__":
    uji_kasus()
    t, sinyal, besar, fs, n_sinyal = studi_kasus()
    hasil = eksperimen()
    buat_grafik(hasil, t, sinyal, besar, fs, n_sinyal)
