#ifndef SERIAL_H
#define SERIAL_H

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int csjac2(const float a11, const float a22, const float a21r, const float a21i, float *const cs, float *const snr, float *const sni, float *const l1, float *const l2);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int zsjac2(const double a11, const double a22, const double a21r, const double a21i, double *const cs, double *const snr, double *const sni, double *const l1, double *const l2);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int wsjac2(const long double a11, const long double a22, const long double a21r, const long double a21i, long double *const cs, long double *const snr, long double *const sni, long double *const l1, long double *const l2);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int ssjac2(const float a11, const float a22, const float a21, float *const cs, float *const sn, float *const l1, float *const l2);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int dsjac2(const double a11, const double a22, const double a21, double *const cs, double *const sn, double *const l1, double *const l2);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int xsjac2(const long double a11, const long double a22, const long double a21, long double *const cs, long double *const sn, long double *const l1, long double *const l2);

#endif /* !SERIAL_H */
