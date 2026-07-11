\[Sigma] = 10;



\[Rho] = 28;



\[Beta] = 8/3;



eqptslns =   NSolve[{0 == \[Sigma] (y - x), 0 == x (\[Rho] - z) - y,   0 == x*y - \[Beta]*z}, {x, y, z}];



eqpts = {x, y, z} /. eqptslns;



(*DXYZ*)
h = 0.01;



xt = Table[1, {k, 0.1, 5000}];



yt = Table[1, {k, 1, Length[xt]}];



zt = Table[1, {k, 1, Length[xt]}];



f[{x_, y_, z_}] = {x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x),
    y +(1-E^-h)/1*((x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))(\[Rho]-z)-y), 
   z +((1-E^(-\[Beta]*h)))/\[Beta]*((x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))*(y +(1-E^-h)/1*((x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))(\[Rho]-z)-y))-\[Beta]*z)};
   
For[k = 1, k <= Length[xt] - 1, k++,
 {xt[[k + 1]], yt[[k + 1]], zt[[k + 1]]} = f[{xt[[k]], yt[[k]], zt[[k]]}];
 (* xt[[k+1]]=xt[[k]]+((1-\[ExponentialE]^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(yt[[k]]-xt[[k]]); *)
 (* yt[[k+1]]=yt[[k]]+(1-\[ExponentialE]^-h)/1*(xt[[k+1]](\[Rho]-zt[[k]])-yt[[k]]); *)
 (* zt[[k+1]]=zt[[k]]+((1-\[ExponentialE]^(-\[Beta]*h)))/\[Beta]*(xt[[k+1]]*yt[[k+1]]-\[Beta]*zt[[k]]); *)
 ]

slnt = Table[{xt[[k]], yt[[k]], zt[[k]]}, {k, 1, Length[xt]}];

Show[{ListPointPlot3D[slnt, PlotStyle -> {PointSize[Medium]}, 
   AxesLabel -> {"x", "y", "z"}, PlotRange -> All], 
  ListPointPlot3D[eqpts, 
   PlotStyle -> {Opacity[0.5], Red, PointSize[0.02]}], 
  Graphics3D[{Red, Line[slnt]}]}, ImageSize -> 600]

jmat[{x_, y_, z_}] = 
  FullSimplify[
   Table[{D[f[{x, y, z}], x][[i]], D[f[{x, y, z}], y][[i]], 
     D[f[{x, y, z}], z][[i]]}, {i, 1, 3}]];

evals = Flatten[
  Table[Eigenvalues[jmat[eqpts[[i]]]], {i, 1, Length[eqpts]}]]

enorms = Table[Norm[evals[[i]]], {i, 1, Length[evals]}]

Show[{ParametricPlot[{Cos[\[Tau]], Sin[\[Tau]]}, {\[Tau], 0, 2 \[Pi]},
    PlotRange -> {{-Max[enorms] - 0.1, 
      Max[enorms] + 0.1}, {-Max[enorms] - 0.1, Max[enorms] + 0.1}}], 
  ListPlot[
  Table[ Tooltip[{Re[evals[[i]]], Im[evals[[i]]]}, {Re[evals[[i]]], Im[evals[[i]]]}], {i, 1, Length[evals]} ], 
   PlotStyle -> {Red}]}, ImageSize -> 400]




(*1. Adım boyutuna (step) bağlı dinamik haritanın tanımlanması*)
fStep[{x_,y_,z_},step_]:={x+((1-E^(-\[Sigma]*step))/\[Sigma])*\[Sigma]*(y-x),y+(1-E^-step)/1*((x+((1-E^(-\[Sigma]*step))/\[Sigma])*\[Sigma]*(y-x)) (\[Rho]-z)-y),z+((1-E^(-\[Beta]*step)))/\[Beta]*((x+((1-E^(-\[Sigma]*step))/\[Sigma])*\[Sigma]*(y-x))*(y+(1-E^-step)/1*((x+((1-E^(-\[Sigma]*step))/\[Sigma])*\[Sigma]*(y-x)) (\[Rho]-z)-y))-\[Beta]*z)};

(*2. slnt yörüngesi üzerinden Richardson Ekstrapolasyonu*)
errors=Table[Module[{pt,ptH,pt2H,gradH,grad2H,errX,errY,errZ},pt=slnt[[k]];
(*h ve 2h ile ileri sarma*)ptH=fStep[pt,h];
pt2H=fStep[pt,2*h];
(*Gradyan hesaplaması*)gradH=(ptH-pt)/h;
grad2H=(pt2H-pt)/(2*h);
(*Hata: |Türev(2h)-Türev(h)|*){errX,errY,errZ}=Abs[grad2H-gradH];
{errX,errY,errZ,errX+errY+errZ}],{k,1,Length[slnt]-1}];

(*3. Verileri Ayrıştırma ve Çizim*)
xErr=errors[[All,1]];
yErr=errors[[All,2]];
zErr=errors[[All,3]];
totalErr=errors[[All,4]];

ListLogPlot[{xErr,yErr,zErr,totalErr},Joined->True,PlotRange->All,PlotStyle->{Blue,Orange,Green,Red},PlotLegends->{"x-error","y-error","z-error","Total error"},AxesLabel->{"iter","abs(w - w')"},PlotLabel->"NSFD [x\[RightArrow]y\[RightArrow]z] Richardson Hata Analizi (h="<>ToString[h]<>")",ImageSize->600]



ClearAll["Global`*"];

(*Parametreler*)
\[Sigma]=10;\[Rho]=28;\[Beta]=8/3;
h=0.01;
iters=5000;

(* ===STANDART EULER YÖNTEMİ===*)
eulerStep[{x_,y_,z_}]:={x+h*\[Sigma]*(y-x),y+h*(x*(\[Rho]-z)-y),z+h*(x*y-\[Beta]*z)};

(* ===SENİN NSFD[x->y->z] YÖNTEMİN===*)
nsfdStep[{x_,y_,z_}]:={x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x),y+(1-E^-h)/1*((x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x)) (\[Rho]-z)-y),z+((1-E^(-\[Beta]*h)))/\[Beta]*((x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))*(y+(1-E^-h)/1*((x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x)) (\[Rho]-z)-y))-\[Beta]*z)};

(*Yörüngeleri Oluşturma*)
pt0={1.,1.,1.};
eulerPts=NestList[eulerStep,pt0,iters];
nsfdPts=NestList[nsfdStep,pt0,iters];

(* ===2. YÖNTEM:GÖRECELİ HATA (RELATIVE ERROR) HESABI===*)
relError=Table[Module[{wc,wq,num,den},wc=eulerPts[[i]];
wq=nsfdPts[[i]];
num=Norm[wc-wq];
den=1+Norm[wc]+Norm[wq];
num/den],{i,1,Length[eulerPts]}];

plotError=ListPlot[relError,Joined->True,PlotStyle->Purple,PlotRange->All,AxesLabel->{"iter","Relative Error"},PlotLabel->"Denklem (29) Goreceli Hata: Euler vs NSFD (h="<>ToString[h]<>")",ImageSize->500];

(* ===3. YÖNTEM:GÖRSEL KARŞILAŞTIRMA (ESTETİK VE YÜKSEK ÇÖZÜNÜRLÜK)===*)(* ===3. YÖNTEM:GÖRSEL KARŞILAŞTIRMA (ÇİZGİ vs NOKTA)===*)
plot3D=Show[(*Euler yörüngesi:Şeffaf mavi çizgiler*)
Graphics3D[{Blue,Opacity[0.5],Thickness[0.001],Line[eulerPts]}],(*NSFD yörüngesi:Belirgin kırmızı noktalar*)Graphics3D[{Red,Opacity[0.9],PointSize[0.003],Point[nsfdPts]}],Axes->True,AxesLabel->{"x","y","z"},(*Karmaşayı azaltan ve kanatları net gösteren kamera açısı*)ViewPoint->{-2,-2.5,1.5},Boxed->True,PlotLabel->Style["Mavi Çizgi: Standart Euler | Kırmızı Noktalar: NSFD",14,Bold],(*Yüksek kalite (1080p dengi) render ayarı*)ImageSize->600];

(*Çıktıları Göster*)
Column[{plotError,plot3D}]

(*Çıktıları Göster*)
Column[{plotError,plot3D}]

(*Çıktıları Göster*)
Column[{plotError,plot3D}]



(*DXZY*)
h = 0.04;



xt = Table[1, {k, 1, 500}];



yt = Table[1, {k, 1, Length[xt]}];



zt = Table[1, {k, 1, Length[xt]}];



f[{x_, y_, z_}] = {x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x),
    y +(1-E^-h)/1*((x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))(\[Rho]-( z +((1-E^(-\[Beta]*h)))/\[Beta]*((x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))*y-\[Beta]*z)))-y), 
   z +((1-E^(-\[Beta]*h)))/\[Beta]*((x +((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))*y-\[Beta]*z)};
   
For[k = 1, k <= Length[xt] - 1, k++,
 {xt[[k + 1]], yt[[k + 1]], zt[[k + 1]]} = f[{xt[[k]], yt[[k]], zt[[k]]}];
 (* xt[[k+1]]=xt[[k]]+((1-\[ExponentialE]^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(yt[[k]]-xt[[k]]); *)
 (* yt[[k+1]]=yt[[k]]+(1-\[ExponentialE]^-h)/1*(xt[[k+1]](\[Rho]-zt[[k]])-yt[[k]]); *)
 (* zt[[k+1]]=zt[[k]]+((1-\[ExponentialE]^(-\[Beta]*h)))/\[Beta]*(xt[[k+1]]*yt[[k+1]]-\[Beta]*zt[[k]]); *)
 ]

slnt = Table[{xt[[k]], yt[[k]], zt[[k]]}, {k, 1, Length[xt]}];

Show[{ListPointPlot3D[slnt, PlotStyle -> {PointSize[Medium]}, 
   AxesLabel -> {"x", "y", "z"}, PlotRange -> All], 
  ListPointPlot3D[eqpts, 
   PlotStyle -> {Opacity[0.5], Red, PointSize[0.02]}], 
  Graphics3D[{Red, Line[slnt]}]}, ImageSize -> 600]

jmat[{x_, y_, z_}] = 
  FullSimplify[
   Table[{D[f[{x, y, z}], x][[i]], D[f[{x, y, z}], y][[i]], 
     D[f[{x, y, z}], z][[i]]}, {i, 1, 3}]];

evals = Flatten[
  Table[Eigenvalues[jmat[eqpts[[i]]]], {i, 1, Length[eqpts]}]]

enorms = Table[Norm[evals[[i]]], {i, 1, Length[evals]}]

Show[{ParametricPlot[{Cos[\[Tau]], Sin[\[Tau]]}, {\[Tau], 0, 2 \[Pi]},
    PlotRange -> {{-Max[enorms] - 0.1, 
      Max[enorms] + 0.1}, {-Max[enorms] - 0.1, Max[enorms] + 0.1}}], 
  ListPlot[
  Table[ Tooltip[{Re[evals[[i]]], Im[evals[[i]]]}, {Re[evals[[i]]], Im[evals[[i]]]}], {i, 1, Length[evals]} ], 
   PlotStyle -> {Red}]}, ImageSize -> 400]




(*DYXZ*)
h = 0.04;



xt = Table[1, {k, 1, 500}];



yt = Table[1, {k, 1, Length[xt]}];



zt = Table[1, {k, 1, Length[xt]}];



f[{x_, y_, z_}] = {x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*((y+(1-E^-h)/1(x(\[Rho]-z)-y))-x),  
y+(1-E^-h)/1*(x(\[Rho]-z)-y), 
   z+((1-E^(-\[Beta]*h))/\[Beta])*((x+(1-E^(-\[Sigma]*h))/\[Sigma]/\[Sigma]*((y+(1-E^-h)/1(x(\[Rho]-z)-y))-x))*(y+(1-E^-h)/1(x(\[Rho]-z)-y))-\[Beta]*z)};
   
For[k = 1, k <= Length[xt] - 1, k++,
 {xt[[k + 1]], yt[[k + 1]], zt[[k + 1]]} = f[{xt[[k]], yt[[k]], zt[[k]]}];
 (* xt[[k+1]]=xt[[k]]+((1-\[ExponentialE]^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(yt[[k]]-xt[[k]]); *)
 (* yt[[k+1]]=yt[[k]]+(1-\[ExponentialE]^-h)/1*(xt[[k+1]](\[Rho]-zt[[k]])-yt[[k]]); *)
 (* zt[[k+1]]=zt[[k]]+((1-\[ExponentialE]^(-\[Beta]*h)))/\[Beta]*(xt[[k+1]]*yt[[k+1]]-\[Beta]*zt[[k]]); *)
 ]

slnt = Table[{xt[[k]], yt[[k]], zt[[k]]}, {k, 1, Length[xt]}];

Show[{ListPointPlot3D[slnt, PlotStyle -> {PointSize[Medium]}, 
   AxesLabel -> {"x", "y", "z"}, PlotRange -> All], 
  ListPointPlot3D[eqpts, 
   PlotStyle -> {Opacity[0.5], Red, PointSize[0.02]}], 
  Graphics3D[{Red, Line[slnt]}]}, ImageSize -> 600]

jmat[{x_, y_, z_}] = 
  FullSimplify[
   Table[{D[f[{x, y, z}], x][[i]], D[f[{x, y, z}], y][[i]], 
     D[f[{x, y, z}], z][[i]]}, {i, 1, 3}]];

evals = Flatten[
  Table[Eigenvalues[jmat[eqpts[[i]]]], {i, 1, Length[eqpts]}]]

enorms = Table[Norm[evals[[i]]], {i, 1, Length[evals]}]

Show[{ParametricPlot[{Cos[\[Tau]], Sin[\[Tau]]}, {\[Tau], 0, 2 \[Pi]},
    PlotRange -> {{-Max[enorms] - 0.1, 
      Max[enorms] + 0.1}, {-Max[enorms] - 0.1, Max[enorms] + 0.1}}], 
  ListPlot[
  Table[ Tooltip[{Re[evals[[i]]], Im[evals[[i]]]}, {Re[evals[[i]]], Im[evals[[i]]]}], {i, 1, Length[evals]} ], 
   PlotStyle -> {Red}]}, ImageSize -> 400]



(*DYZX*)
h = 0.04;



xt = Table[1, {k, 1, 500}];



yt = Table[1, {k, 1, Length[xt]}];



zt = Table[1, {k, 1, Length[xt]}];



f[{x_, y_, z_}] = {x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*((y+(1-E^-h)/1(x(\[Rho]-z)-y))-x),  
y+(1-E^-h)/1*(x(\[Rho]-z)-y), 
   z+((1-E^(-\[Beta]*h))/\[Beta])*(x*(y+(1-E^-h)/1(x(\[Rho]-z)-y))-\[Beta]*z)};
   
For[k = 1, k <= Length[xt] - 1, k++,
 {xt[[k + 1]], yt[[k + 1]], zt[[k + 1]]} = f[{xt[[k]], yt[[k]], zt[[k]]}];
 (* xt[[k+1]]=xt[[k]]+((1-\[ExponentialE]^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(yt[[k]]-xt[[k]]); *)
 (* yt[[k+1]]=yt[[k]]+(1-\[ExponentialE]^-h)/1*(xt[[k+1]](\[Rho]-zt[[k]])-yt[[k]]); *)
 (* zt[[k+1]]=zt[[k]]+((1-\[ExponentialE]^(-\[Beta]*h)))/\[Beta]*(xt[[k+1]]*yt[[k+1]]-\[Beta]*zt[[k]]); *)
 ]

slnt = Table[{xt[[k]], yt[[k]], zt[[k]]}, {k, 1, Length[xt]}];

Show[{ListPointPlot3D[slnt, PlotStyle -> {PointSize[Medium]}, 
   AxesLabel -> {"x", "y", "z"}, PlotRange -> All], 
  ListPointPlot3D[eqpts, 
   PlotStyle -> {Opacity[0.5], Red, PointSize[0.02]}], 
  Graphics3D[{Red, Line[slnt]}]}, ImageSize -> 600]

jmat[{x_, y_, z_}] = 
  FullSimplify[
   Table[{D[f[{x, y, z}], x][[i]], D[f[{x, y, z}], y][[i]], 
     D[f[{x, y, z}], z][[i]]}, {i, 1, 3}]];

evals = Flatten[
  Table[Eigenvalues[jmat[eqpts[[i]]]], {i, 1, Length[eqpts]}]]

enorms = Table[Norm[evals[[i]]], {i, 1, Length[evals]}]

Show[{ParametricPlot[{Cos[\[Tau]], Sin[\[Tau]]}, {\[Tau], 0, 2 \[Pi]},
    PlotRange -> {{-Max[enorms] - 0.1, 
      Max[enorms] + 0.1}, {-Max[enorms] - 0.1, Max[enorms] + 0.1}}], 
  ListPlot[
  Table[ Tooltip[{Re[evals[[i]]], Im[evals[[i]]]}, {Re[evals[[i]]], Im[evals[[i]]]}], {i, 1, Length[evals]} ], 
   PlotStyle -> {Red}]}, ImageSize -> 400]



(*DZXY*)
h = 0.04;



xt = Table[1, {k, 1, 500}];



yt = Table[1, {k, 1, Length[xt]}];



zt = Table[1, {k, 1, Length[xt]}];



f[{x_, y_, z_}] = {x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x),  
y+(1-E^-h)/1*((x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(y-x))(\[Rho]-(z+((1-E^(-\[Beta]*h))/\[Beta])*(x*y-\[Beta]*z)))-y), 
   z+((1-E^(-\[Beta]*h))/\[Beta])*(x*y-\[Beta]*z)};
   
For[k = 1, k <= Length[xt] - 1, k++,
 {xt[[k + 1]], yt[[k + 1]], zt[[k + 1]]} = f[{xt[[k]], yt[[k]], zt[[k]]}];
 (* xt[[k+1]]=xt[[k]]+((1-\[ExponentialE]^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(yt[[k]]-xt[[k]]); *)
 (* yt[[k+1]]=yt[[k]]+(1-\[ExponentialE]^-h)/1*(xt[[k+1]](\[Rho]-zt[[k]])-yt[[k]]); *)
 (* zt[[k+1]]=zt[[k]]+((1-\[ExponentialE]^(-\[Beta]*h)))/\[Beta]*(xt[[k+1]]*yt[[k+1]]-\[Beta]*zt[[k]]); *)
 ]

slnt = Table[{xt[[k]], yt[[k]], zt[[k]]}, {k, 1, Length[xt]}];

Show[{ListPointPlot3D[slnt, PlotStyle -> {PointSize[Medium]}, 
   AxesLabel -> {"x", "y", "z"}, PlotRange -> All], 
  ListPointPlot3D[eqpts, 
   PlotStyle -> {Opacity[0.5], Red, PointSize[0.02]}], 
  Graphics3D[{Red, Line[slnt]}]}, ImageSize -> 600]

jmat[{x_, y_, z_}] = 
  FullSimplify[
   Table[{D[f[{x, y, z}], x][[i]], D[f[{x, y, z}], y][[i]], 
     D[f[{x, y, z}], z][[i]]}, {i, 1, 3}]];

evals = Flatten[
  Table[Eigenvalues[jmat[eqpts[[i]]]], {i, 1, Length[eqpts]}]]

enorms = Table[Norm[evals[[i]]], {i, 1, Length[evals]}]

Show[{ParametricPlot[{Cos[\[Tau]], Sin[\[Tau]]}, {\[Tau], 0, 2 \[Pi]},
    PlotRange -> {{-Max[enorms] - 0.1, 
      Max[enorms] + 0.1}, {-Max[enorms] - 0.1, Max[enorms] + 0.1}}], 
  ListPlot[
  Table[ Tooltip[{Re[evals[[i]]], Im[evals[[i]]]}, {Re[evals[[i]]], Im[evals[[i]]]}], {i, 1, Length[evals]} ], 
   PlotStyle -> {Red}]}, ImageSize -> 400]



(*DZYX*)
h = 0.04;



xt = Table[1, {k, 1, 500}];



yt = Table[1, {k, 1, Length[xt]}];



zt = Table[1, {k, 1, Length[xt]}];



f[{x_, y_, z_}] = {x+((1-E^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*((y+(1-E^-h)/1*(x(\[Rho]-(z+((1-E^(-\[Beta]*h))/\[Beta])*(x*y-\[Beta]*z)))-y))-x),  
y+(1-E^-h)/1*(x(\[Rho]-(z+((1-E^(-\[Beta]*h))/\[Beta])*(x*y-\[Beta]*z)))-y), 
   z+((1-E^(-\[Beta]*h))/\[Beta])*(x*y-\[Beta]*z)};
   
For[k = 1, k <= Length[xt] - 1, k++,
 {xt[[k + 1]], yt[[k + 1]], zt[[k + 1]]} = f[{xt[[k]], yt[[k]], zt[[k]]}];
 (* xt[[k+1]]=xt[[k]]+((1-\[ExponentialE]^(-\[Sigma]*h))/\[Sigma])*\[Sigma]*(yt[[k]]-xt[[k]]); *)
 (* yt[[k+1]]=yt[[k]]+(1-\[ExponentialE]^-h)/1*(xt[[k+1]](\[Rho]-zt[[k]])-yt[[k]]); *)
 (* zt[[k+1]]=zt[[k]]+((1-\[ExponentialE]^(-\[Beta]*h)))/\[Beta]*(xt[[k+1]]*yt[[k+1]]-\[Beta]*zt[[k]]); *)
 ]

slnt = Table[{xt[[k]], yt[[k]], zt[[k]]}, {k, 1, Length[xt]}];

Show[{ListPointPlot3D[slnt, PlotStyle -> {PointSize[Medium]}, 
   AxesLabel -> {"x", "y", "z"}, PlotRange -> All], 
  ListPointPlot3D[eqpts, 
   PlotStyle -> {Opacity[0.5], Red, PointSize[0.02]}], 
  Graphics3D[{Red, Line[slnt]}]}, ImageSize -> 600]

jmat[{x_, y_, z_}] = 
  FullSimplify[
   Table[{D[f[{x, y, z}], x][[i]], D[f[{x, y, z}], y][[i]], 
     D[f[{x, y, z}], z][[i]]}, {i, 1, 3}]];

evals = Flatten[
  Table[Eigenvalues[jmat[eqpts[[i]]]], {i, 1, Length[eqpts]}]]

enorms = Table[Norm[evals[[i]]], {i, 1, Length[evals]}]

Show[{ParametricPlot[{Cos[\[Tau]], Sin[\[Tau]]}, {\[Tau], 0, 2 \[Pi]},
    PlotRange -> {{-Max[enorms] - 0.1, 
      Max[enorms] + 0.1}, {-Max[enorms] - 0.1, Max[enorms] + 0.1}}], 
  ListPlot[
  Table[ Tooltip[{Re[evals[[i]]], Im[evals[[i]]]}, {Re[evals[[i]]], Im[evals[[i]]]}], {i, 1, Length[evals]} ], 
   PlotStyle -> {Red}]}, ImageSize -> 400]
