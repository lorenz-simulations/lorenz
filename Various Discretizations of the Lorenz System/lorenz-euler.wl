(*V1*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {x + h*\[Sigma] (y - x), y + h (x (\[Rho] - z) - y), z + h (x*y - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]

(*V2*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {(x + h*\[Sigma]*y)/(1 + h*\[Sigma]),y + h*((x + h*\[Sigma]*y)/(1 + h*\[Sigma]) (\[Rho] - z) - y),z + h*((x + h*\[Sigma]*y)/(1 + h*\[Sigma]) y - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]

(*V3*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {x +  h*\[Sigma] ((y + h*x (\[Rho] - z))/(1 + h) - x), ( y + h*x (\[Rho] - z))/(1 + h), z + h (x*(y + h*x (\[Rho] - z))/(1 + h) - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts],Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]

(*V4*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {x + h*\[Sigma] (y - x), y + h (x (\[Rho] -( z + h x*y )/(1+ \[Beta]*h)) - y), ( z + h x*y )/(1+ \[Beta]*h)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]
 
(*V2&3*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {(x*(1+h)+h*\[Sigma] *y)/((1+h*\[Sigma] )(1+h)-h^2*(\[Sigma] *(\[Rho]-z))), 
(y*(1+h*\[Sigma])+h*x *(\[Rho]-z))/((1+h*\[Sigma] )(1+h)-h^2*(\[Sigma] *(\[Rho]-z))), 
z + h *((x*(1+h)+h*\[Sigma] *y)/((1+h*\[Sigma] )*(1+h)-h^2*(\[Sigma] *(\[Rho]-z)))*(y*(1+h*\[Sigma])+h*x *(\[Rho]-z))/((1+h*\[Sigma] )*(1+h)-h^2(\[Sigma] *(\[Rho]-z)))- \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]
 
(*V2&4*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {(x + h*\[Sigma]*y)/(1 + h*\[Sigma]),
y + h ((x + h*\[Sigma]*y)/(1 + h*\[Sigma]) (\[Rho] - ( z + h x*y )/(1+ \[Beta]*h)) - y),
( z + h x*y )/(1+ \[Beta]*h)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]

(*V3&4*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := 
{x + h*\[Sigma] ((z*(1+h*x)-(1+h*\[Beta])*(y+h*\[Rho]*x))/((1+h)*(1+h*\[Beta])+h*x*(1+h*x)) - x),
(z*(1+h*x)-(1+h*\[Beta])*(y+h*\[Rho]*x))/((1+h)*(1+h*\[Beta])+h*x*(1+h*x)), 
-(z+h*(z+x*(y+h*\[Rho]*x)))/((1+h)*(1+h*\[Beta])+h*x*(1+h*x))}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]

(*V5*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {x + h*\[Sigma] (y - x),
y + h ((x + h*\[Sigma]*y)/(1 + h*\[Sigma])* (\[Rho] - z) - y), 
z + h ((x + h*\[Sigma]*y)/(1 + h*\[Sigma])* ( y + h*x (\[Rho] - z))/(1 + h) - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]
 
(*V6*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := {x + h*\[Sigma] (y - x),
y + h ((x + h*\[Sigma]*y)/(1 + h*\[Sigma]) (\[Rho] -  ( z + h x*y )/(1+ \[Beta]*h)) - y), 
z + h ((x + h*\[Sigma]*y)/(1 + h*\[Sigma])*y - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]
 
(*V7*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01; 
f[{x_, y_, z_}] := 
{ x + h*\[Sigma]*((y + h*(x*(\[Rho] - z) - y)) - x),
y + h*(x*(\[Rho] - z) - y),
z + h*((x + h*\[Sigma]*((y + h*(x*(\[Rho] - z) - y)) - x))*(y + h*(x*(\[Rho] - z) - y)) - \[Beta]*z) }
pts = NestList[f, {1, 1, 1}, 5000]; 
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]

(*V8*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := 
{x + h*\[Sigma]*((y + h*(x*(\[Rho] - z) - y)) - x),
  y + h*(x*(\[Rho] - z) - y),
  z + h*(x*(y + h*(x*(\[Rho] - z) - y)) - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]
 
(*V9*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] := 
{x + h*\[Sigma]*(y - x),
  y + h*((x + h*\[Sigma]*(y - x))*(\[Rho] - (z + h*(x*y - \[Beta]*z))) - y),
  z + h*(x*y - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]
 
(*V10*)
\[Sigma] = 10;
\[Rho] = 28;
\[Beta] = 8/3;
h = 0.01;
f[{x_, y_, z_}] :=
{x + h*\[Sigma]*((y + h*(x*(\[Rho] - (z + h*(x*y - \[Beta]*z))) - y)) - x),
  y + h*(x*(\[Rho] - (z + h*(x*y - \[Beta]*z))) - y),
  z + h*(x*y - \[Beta]*z)}
pts = NestList[f, {1, 1, 1}, 5000];
Show[{ListPointPlot3D[pts], Graphics3D[{Opacity[0.35], Red, Line[pts]}]}, ViewPoint -> Above]
