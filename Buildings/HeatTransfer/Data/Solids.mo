within Buildings.HeatTransfer.Data;
package Solids
  "Package with solid material, characterized by thermal conductance, density and specific heat capacity"
  extends Modelica.Icons.MaterialPropertiesPackage;

  record Generic "Thermal properties of solids with heat storage"
    extends Buildings.HeatTransfer.Data.BaseClasses.Material(
      final R=x/k,
      final TSol=293.15,
      final TLiq=293.15,
      final LHea=0,
      final phasechange=false);
    annotation (defaultComponentName="mat", Documentation(info="<html>
<p>
Generic record for solid materials.
The material is characterized by its 
thermal conductivity, mass density and specific
heat capacity.
</p>
</html>", revisions="<html>
<ul>
<li>
September 9, 2010, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
  end Generic;

  record Brick = Buildings.HeatTransfer.Data.Solids.Generic (
      k=0.89,
      d=1920,
      c=790) "Brick (k=0.89)";
  record Concrete = Buildings.HeatTransfer.Data.Solids.Generic (
      k=1.4,
      d=2240,
      c=840) "Concrete (k=1.4)";
  record InsulationBoard = Buildings.HeatTransfer.Data.Solids.Generic (
      k=0.03,
      d=40,
      c=1200) "Insulation board (k=0.03)";
  record GypsumBoard = Buildings.HeatTransfer.Data.Solids.Generic (
      k=0.16,
      d=800,
      c=1090) "Gypsum board (k=0.58)";
  record Plywood = Buildings.HeatTransfer.Data.Solids.Generic (
      k=0.12,
      d=540,
      c=1210) "Plywood (k=0.12)";
  record GenericHM
    extends BaseClasses.HygroThermalMaterial(
      final R=x/k,
      final TSol=293.15,
      final TLiq=293.15,
      final LHea=0,
      final phasechange=false);
    annotation (defaultComponentName="mat", Documentation(info="<html>
<p>
Generic record for solid materials.
The material is characterized by its 
thermal conductivity, mass density and specific
heat capacity.
</p>
</html>", revisions="<html>
<ul>
<li>
September 9, 2010, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));

  end GenericHM;

  record ConcreteHM = Buildings.HeatTransfer.Data.Solids.GenericHM (
      k=1.6,
      d=2300,
      c=850,
      w_80=85,
      w_f=150,
      mu=6,
      A_layer=0.003,
      sorp_tab={{0.0,0.0},{0.5,48.0},{0.6,58.0},{0.8,85.0},{0.9,100.0},{0.95,
          118.0},{1.0,150.0}},
      dws_tab=[0.0, 1.0e-16; 72.0, 7.4e-11; 85.0, 2.5e-10; 100.0, 1.0e-9; 118.0,
          1.2e-9],
      dww_tab=[0.0, 1.0e-16; 72.0, 7.4e-12; 85.0, 2.5e-11; 100.0, 1.0e-10;
          118.0, 1.2e-10],
      my_tab=[0.0, 180; 1.0, 180],
      lamb_tab=[0.0, 1.6; 180, 1.6],
      b_h=8,
      por=0.18);
  record GypsumBoardHM = Buildings.HeatTransfer.Data.Solids.GenericHM (
      k=0.2,
      d=850,
      c=850,
      w_80=6.3,
      w_f=400,
      mu=9,
      por=0.65,
      b_h=4,
      A_layer=0.287,
      sorp_tab=[0.0, 0.0; 0.5, 3.6; 0.65, 5.2; 0.8, 6.3; 0.9, 11.0; 0.93, 17.0;
          0.95, 19.0; 0.99, 113.0; 0.995, 124.0; 0.999, 328.0; 0.9995, 378.0;
          1.0, 400.0],
      lamb_tab=[0.0, 0.2; 650.0, 0.2],
      dws_tab=[0.0, 1.0e-16; 60, 3.0e-9; 100, 1.0e-7; 160, 1.0e-7; 240, 1.2e-7;
          320, 2.2e-7; 360, 6.0e-7; 380, 9.0e-7; 400, 4.5e-6],
      dww_tab=[0.0, 1.0e-16; 60, 3.0e-9; 100, 8.0e-9; 160, 8.0e-9; 240, 1.3e-8;
          320, 1.0e-7; 360, 3.0e-7; 380, 7.0e-7; 400, 1.0e-6],
      my_tab=[0.0, 9.0; 1.0, 9.0]);
  annotation (Documentation(info="<html>
<p>
Package with records for solid materials.
The material is characterized by its 
thermal conductivity, mass density and specific
heat capacity.
</p>
<p>
These material records automatically compute the spatial grid
that is used to compute transient heat conduction.
In building materials, the thermal diffusivity of adjacent layer materials can differ by an order of magnitude. If the spatial grid generation were not to account for the material properties, then the time rate of change of the different temperature nodes would be significantly different from each other.
Therefore, records in the packages
<a href=\"Buildings.HeatTransfer.Data.Solids\">
Buildings.HeatTransfer.Data.Solids</a>
and
<a href=\"Buildings.HeatTransfer.Data.SolidsPCM\">
Buildings.HeatTransfer.Data.SolidsPCM</a>
generate the spatial grid so that under the assumption of equal heat transfer, each node temperature has a similar time rate of change.
</p>
<p>
The computation is as follows:
</p>
<p>
From dimensionless analysis, one can obtain a characteristic time, called the <em>Fourier</em> number, as
</p>
<p align=\"center\" style=\"font-style:italic;\">
Fo = &alpha; t &frasl; L<sup>2</sup>
</p>
<p>
where <i>&alpha;</i> denotes the thermal diffusivity, <i>t</i> denotes time and <i>L</i> denotes the characteristic length. 
We like to generate the spatial grid so that the ratio
<i>t &frasl; Fo</i>
is equal to an arbitrary constant 
<i>&Pi;</i>, which we define as
</p>
<p align=\"center\" style=\"font-style:italic;\">
&Pi; = ( t &frasl; Fo )<sup>1/2</sup> 
</p>

<p>and hence</p>

<p align=\"center\" style=\"font-style:italic;\">
&Pi; = L &frasl; &radic; &alpha;.
</p>

<p>
Now, let <i>x</i>
denote the thickness of the material layer.
Then, we compute the time constant of the material layer as
</p>
<p align=\"center\" style=\"font-style:italic;\">
&Pi;<sub>x</sub> = x &frasl; &radic; &alpha;,
</p>
<p>
and we compute the estimated number of elements <i>N' &isin; &#8477;</i> 
for the material layer as</p>

<p align=\"center\" style=\"font-style:italic;\">
N' = N<sub>ref</sub> &Pi;<sub>x</sub> &frasl; &Pi;<sub>ref</sub>
</p>

<p>
where <i>&Pi;<sub>ref</sub> &isin; &#8469;</i> is a user-specified number of elements 
for a reference material, which is equal to the parameter
<code>nStaRef</code>, and defined as a concrete construction with thickness 
<i>L<sub>ref</sub> = 0.20</i> meter and thermal diffusivity
<i>&alpha;<sub>ref</sub> = 3.64E-7</i> m<sup>2</sup>/s.
Hence, 
<i>&Pi;<sub>ref</sub> = L<sub>ref</sub>/ &radic; &alpha;<sub>ref</sub> = 331.4</i>
&radic;s.
</p>

<p>
Next, we define the number of elements for the material layer as
<p align=\"center\" style=\"font-style:italic;\">
<i>N<sub>x</sub> = &lceil;  N' &rceil;</i>
</p>

<p>
where the notation <i>&lceil; &#8901; &rceil;</i> is defined, for
<i>s &isin; &#8477;</i>, as
<p align=\"center\" style=\"font-style:italic;\">
&lceil; s &rceil; = min{ k &isin; &#8484; | k &ge; s }.
</p>
<p>
Finally, we divide the material layer in compartments of length 
<i>&Delta; = x &frasl; N<sub>x</sub></i>.
</p>

</html>", revisions="<html>
<ul>
<li>
September 9, 2010, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
end Solids;
