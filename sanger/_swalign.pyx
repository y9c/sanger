# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: language = c
"""
Cython Smith-Waterman local aligner (self-contained; no external ssw, no NumPy).

Compiled into ``sanger._swalign`` when a C compiler is available at build time;
otherwise sanger falls back to the NumPy-vectorised / pure-Python
implementation in :mod:`sanger.align`.

``sw_align(reference, query)`` returns ``(reference_begin, query_begin,
query_aligned, ref_aligned)`` where the aligned strings use ``-`` for gaps,
matching the interface consumed by :func:`sanger.align.run_align`.

The DP matrices are allocated with ``malloc`` so the extension needs only the C
compiler -- no NumPy headers -- keeping builds light and portable.
"""

from libc.stdlib cimport malloc, free

cdef int _align(const unsigned char* ref, int n,
                const unsigned char* qry, int m,
                int match, int mismatch, int gap,
                int* dp, unsigned char* tr,
                int* out_bi, int* out_bj, int* out_best) nogil:
    cdef int i, j, up, diag, left, best_local, move, sc, stride = m + 1
    cdef int best = 0, bi = 0, bj = 0
    cdef unsigned char rc
    for i in range(1, n + 1):
        rc = ref[i - 1]
        for j in range(1, m + 1):
            sc = match if rc == qry[j - 1] else mismatch
            diag = dp[(i - 1) * stride + j - 1] + sc
            up = dp[(i - 1) * stride + j] + gap
            left = dp[i * stride + j - 1] + gap
            best_local = diag
            move = 1
            if up > best_local:
                best_local = up; move = 2
            if left > best_local:
                best_local = left; move = 3
            if best_local < 0:
                best_local = 0; move = 0
            dp[i * stride + j] = best_local
            tr[i * stride + j] = move
            if best_local > best:
                best = best_local; bi = i; bj = j
    out_bi[0] = bi
    out_bj[0] = bj
    out_best[0] = best
    return 0


def sw_align(str reference, str query, int match=2, int mismatch=-1, int gap=-1):
    """Local Smith-Waterman alignment of ``query`` onto ``reference``."""
    ref_b = reference.encode("ascii")
    qry_b = query.encode("ascii")
    cdef int n = len(ref_b), m = len(qry_b)
    cdef int stride = m + 1
    if n == 0 or m == 0:
        return 0, 0, query, reference
    cdef const unsigned char* refp = ref_b
    cdef const unsigned char* qryp = qry_b
    cdef int* dp = <int*> malloc((n + 1) * stride * sizeof(int))
    cdef unsigned char* tr = <unsigned char*> malloc((n + 1) * stride * sizeof(unsigned char))
    if dp == NULL or tr == NULL:
        if dp != NULL: free(dp)
        if tr != NULL: free(tr)
        raise MemoryError()
    cdef int idx, total = (n + 1) * stride
    for idx in range(total):
        dp[idx] = 0
        tr[idx] = 0
    cdef int bi = 0, bj = 0, best = 0
    with nogil:
        _align(refp, n, qryp, m, match, mismatch, gap, dp, tr, &bi, &bj, &best)

    # traceback
    cdef int i = bi, j = bj, move
    qal = []
    ral = []
    while i > 0 and j > 0:
        move = tr[i * stride + j]
        if move == 0:
            break
        if move == 1:
            qal.append(chr(qry_b[j - 1]))
            ral.append(chr(ref_b[i - 1]))
            i -= 1; j -= 1
        elif move == 2:
            qal.append("-")
            ral.append(chr(ref_b[i - 1]))
            i -= 1
        else:
            qal.append(chr(qry_b[j - 1]))
            ral.append("-")
            j -= 1
    free(dp)
    free(tr)
    return i, j, "".join(reversed(qal)), "".join(reversed(ral))
