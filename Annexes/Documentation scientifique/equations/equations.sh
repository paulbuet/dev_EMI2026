#!/bin/bash

# Generate equations in png 1/2

tex2png.py '$N= C \lambda^x$' N_1M.png
tex2png.py '$m(D) = aD^b$' m.png
tex2png.py '$v(D) = cD^d \, (\rho_{00}/\rho_{dref})^{0.4}$' v.png
tex2png.py '$n(D)dD = N\,g(D)dD$' n.png
tex2png.py '$g(D)=\frac{\displaystyle{\alpha}}{\displaystyle{\Gamma(\nu)}}\lambda^{\alpha \nu} D ^{\alpha \nu -1} \exp\big(-(\lambda D)^{\alpha}\big)$' g.png
tex2png.py '$M(p)=\int^{\infty}_{0} \, D^{p} g(D) \, dD=\frac{\displaystyle{G(p)}}{\displaystyle{\lambda^{p}}}$' M.png
tex2png.py '$G(p) = \frac{\displaystyle{\Gamma(\nu+p/\alpha)}}{\displaystyle{\Gamma(\nu)}}$' G.png
tex2png.py '$\rho r=\int^{\infty}_{0} \, m(D) n(D) \, dD=\int^{\infty}_{0} \, aD^b N g(D) dD$' rho_r_1.png
tex2png.py '$\rho r=a N M(b)=aN\frac{G(b)}{\lambda^b}$' rho_r_2.png
tex2png.py '$\rho r=aCG(b) \lambda^{x-b}$' rho_r_3_1M.png
tex2png.py '$\rho r=aNG(b) \lambda^{-b}$' rho_r_3_2M.png
tex2png.py '$\lambda = \Big(\frac{\displaystyle{\rho r}}{\displaystyle{aCG(b)}}\Big)^{\frac{\displaystyle{1}}{\displaystyle{x-b}}}$' lambda_1M.png
tex2png.py '$\lambda = \Big(\frac{\displaystyle{\rho r}}{\displaystyle{aNG(b)}}\Big)^{\frac{\displaystyle{1}}{\displaystyle{-b}}}$' lambda_2M.png

# Generate equations in png 2/2

tex2png.py '$\frac{d\rho r}{dt}=\frac{\partial F}{\partial z}$' sedim1.png
tex2png.py '$F=\int^{\infty}_{0} \, v(D)m(D)n(D)dD$' sedim2.png
tex2png.py '$F=\int^{\infty}_{0} \, cD^d \, (\frac{\rho_{00}}{\rho_{dref}})^{0.4}     aD^b    N\,g(D)dD$' sedim3.png
tex2png.py '$F=Nac (\frac{\rho_{00}}{\rho_{dref}})^{0.4} \frac{G(b+d)}{\lambda^{b+d}}$' sedim4.png
tex2png.py '$F=Cac (\frac{\rho_{00}}{\rho_{dref}})^{0.4} G(b+d) (\frac{\rho r}{aCG(b)})^{\frac{b+d-x}{b-x}}$' sedim5_1M.png
tex2png.py '$F=Nac (\frac{\rho_{00}}{\rho_{dref}})^{0.4} G(b+d) (\frac{\rho r}{aNG(b)})^{\frac{b+d}{b}}$' sedim5_2M.png
tex2png.py '$F=\int^{\infty}_{0} \, v(D)n(D)dD$' sedim6.png

# Generate parameter table

tex2png.py '
%\begin{table}
%\caption{Set of parameters used to characterize each ice category
%and the raindrops (Kessler scheme).}
%\begin{center}\label{table2}
\begin{tabular}{|l|l|l|l|l|l|l|}
\hline
Parameters & $r_i$ & $r_s$ & $r_g$ && $r_r$ & $r_c$ \\
\hline \hline
$\alpha$ & 3 & 1 & 1 && 1 & 3 on sea; 1 on land \\
$\nu$    & 3 & 1 & 1 && 1 & 1 on sea; 3 on land \\
\hline
$a$ & 0.82 & 0.02 & 19.6 && 524 & 524 \\
$b$ & 2.5  & 1.9  & 2.8 && 3 & 3  \\
\hline
$c$ & 800  & 5.1  & 124 && 842 & 3.2 10$^7$ \\
$d$ & 1.00 & 0.27 & 0.66 && 0.8 & 2 \\
\hline
$C$ & & 5 & 5 10$^5$ && 8 10$^6$ & \\
$x$ & & 1 &    -0.5   &&   -1  &  \\

\hline
$\overline{f}_0$ & 1.00 & 0.86 & 0.86 && 1.00 & \\
$\overline{f}_1$ &      & 0.28 & 0.28 && 0.26 & \\
$\overline{f}_2$ & 0.14 &      &      &&      & \\
\hline
${\cal C}_{1}$ & $1/\pi$ & 1/$\pi$ & 0.5 && 0.5 & \\
\hline
\end{tabular}
%\end{center}
%\end{table}
' table.png
