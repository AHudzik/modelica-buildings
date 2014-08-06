within Buildings.HeatTransfer.Convection;
model ExteriorHM "Model for an exterior convective heat and mass transfer"
  extends BaseClasses.PartialConvectionHM;

  parameter Buildings.HeatTransfer.Types.ExteriorConvection conMod=
    Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind
    "Convective heat transfer model"
  annotation(Evaluate=true);
  parameter Buildings.HeatTransfer.Types.SurfaceRoughness roughness=
    Buildings.HeatTransfer.Types.SurfaceRoughness.Medium "Surface roughness"
    annotation (Dialog(enable=(conMod <> Buildings.HeatTransfer.Types.InteriorConvection.Fixed)));

  parameter Real beta = 75E-9 "Indoor water vapour coefficient (Kunzel)";
  constant Real hv = 2500000 "Latent heat of phase change";
  parameter Modelica.SIunits.Angle azi "Surface azimuth";

Modelica.Blocks.Interfaces.RealInput v(unit="m/s") "Wind speed"
    annotation (Placement(transformation(extent={{-140,80},{-100,120}})));
  Modelica.Blocks.Interfaces.RealInput dir(unit="rad", displayUnit="deg",
     min=0, max=2*Modelica.Constants.pi) "Wind direction (0=wind from North)"
    annotation (Placement(transformation(extent={{-140,30},{-100,70}})));
  Modelica.SIunits.CoefficientOfHeatTransfer hF
    "Convective heat transfer coefficient due to forced convection";
  Modelica.SIunits.HeatFlux qN_flow
    "Convective heat flux from solid -> fluid due to natural convection";
  Modelica.SIunits.HeatFlux qF_flow
    "Convective heat flux from solid -> fluid due to forced convection";

protected
   parameter Real R(fixed=false) "Surface roughness";
   Real W(min=0.5, max=1) "Wind direction modifier";
initial equation
  if (roughness == Buildings.HeatTransfer.Types.SurfaceRoughness.VeryRough) then
    R=2.17;
  elseif (roughness == Buildings.HeatTransfer.Types.SurfaceRoughness.Rough) then
    R=1.67;
  elseif (roughness == Buildings.HeatTransfer.Types.SurfaceRoughness.Medium) then
    R=1.52;
  elseif (roughness == Buildings.HeatTransfer.Types.SurfaceRoughness.MediumSmooth) then
    R=1.13;
  elseif (roughness == Buildings.HeatTransfer.Types.SurfaceRoughness.Smooth) then
    R=1.11;
  elseif (roughness == Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth) then
    R=1.00;
  else
    R=0;
  end if;

equation
  if (conMod == Buildings.HeatTransfer.Types.ExteriorConvection.Fixed) then
    qN_flow = hFixed * dT;
    W = 1;
    hF = 0;
    qF_flow = 0;
  else
    // Even if hCon is a step function with a step at zero,
    // the product hCon*dT is differentiable at zero with
    // a continuous first derivative
    if isCeiling then
       qN_flow = Buildings.HeatTransfer.Convection.Functions.HeatFlux.ceiling(
                                                                             dT=dT);
    elseif isFloor then
       qN_flow = Buildings.HeatTransfer.Convection.Functions.HeatFlux.floor(
                                                                           dT=dT);
    else
       qN_flow = Buildings.HeatTransfer.Convection.Functions.HeatFlux.wall(
                                                                          dT=dT);
    end if;
    // Forced convection
    W = Buildings.HeatTransfer.Convection.Functions.windDirectionModifier(
                                                               azi=azi, dir=dir);
    hF = 2.537 * W * R * 2 / A^(0.25) *
         Buildings.Utilities.Math.Functions.regNonZeroPower(x=v, n=0.5, delta=0.5);
    qF_flow = hF*dT;
  end if;
  q_flow = qN_flow + qF_flow;
  g_flow = beta * dP + m_flow *hv;
  annotation (Documentation(info="<html>
<p>
This is a model for a convective heat and mass transfer for exterior, outside-facing surfaces.
The parameter <code>conMod</code> determines the model that is used to compute
the heat transfer coefficient:
</p>

<ol>
<li><p>If <code>conMod=
<a href=\"modelica://Buildings.HeatTransfer.Types.ExteriorConvection\">
Buildings.HeatTransfer.Types.ExteriorConvection.Fixed</a>
</code>, then
the convective heat transfer coefficient is set to the value specified by the parameter
<code>hFixed</code>.
</p>
</li>
<li>
<p>
If <code>conMod=
<a href=\"modelica://Buildings.HeatTransfer.Types.ExteriorConvection\">
Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind</a>
</code>,
then the convective heat transfer coefficient is
computed based on wind speed, wind direction and temperature difference.
</p>
<p>
The total convection coefficient <i>h<sub>t</sub></i> is the sum of the 
temperature-driven free convection coefficient <i>h<sub>n</sub></i>
and the wind-driven forced convection coefficient <i>h<sub>f</sub></i>,
<p align=\"center\" style=\"font-style:italic;\">
 h<sub>t</sub> = h<sub>n</sub> + h<sub>f</sub>
</p>
The free convection coefficient <i>h<sub>n</sub></i> is computed in the same way as in 
<a href=\"modelica://Buildings.HeatTransfer.Convection.Interior\">
Buildings.HeatTransfer.Convection.Interior</a>.
The forced convection coefficient <i>h<sub>f</sub></i>
is computed based on a correlation by Sparrow, Ramsey, and Mass
(1979), which is 
<p align=\"center\" style=\"font-style:italic;\">
 h<sub>f</sub> = 2.537 W R &radic;( P v &frasl; A )
</p>
<p>
where <i>W=1</i> for windward surfaces and 
<i>W=0.5</i> for leeward surfaces, with leeward defined as greater than 100 degrees
from normal incidence,
<i>R</i> is a surface roughness multiplier,
<i>P</i> is the perimeter of the surface and
<i>A</i> is the area of the surface.
This is the same equation as implemented in EnergyPlus 6.0.
</p>
<p>
We make the simplified assumption that the surface is square, and hence we set
<p align=\"center\" style=\"font-style:italic;\">
 h<sub>f</sub> = 2.537 W R &radic;( 4 v &frasl; &radic;(A) )
</p>
<p>
The surface roughness is specified by the parameter <code>surfaceRoughness</code>
which has to be set to a type of
<a href=\"modelica://Buildings.HeatTransfer.Types.SurfaceRoughness\">
Buildings.HeatTransfer.Types.SurfaceRoughness</a>.The coefficients for the surface roughness are
</p>

<table summary=\"summary\" border=\"1\">
<tr>
<th>Roughness index</th>
<th><i>R</i></th>
<th>Example material</th>
</tr>
<tr><td>VeryRough</td>   <td>2.17</td>  <td>Stucco</td></tr>
<tr><td>Rough</td>        <td>1.67</td>  <td>Brick</td></tr>
<tr><td>MediumRough</td> <td>1.52</td>  <td>Concrete</td></tr>
<tr><td>MediumSmooth</td><td>1.13</td>  <td>Clear pine</td></tr>
<tr><td>Smooth</td>       <td>1.11</td>  <td>Smooth plaster</td></tr>
<tr><td>VerySmooth</td>  <td>1.00</td>  <td>Glass</td></tr>
</table>
The water vapour transfer is calculated in a manner similar to the heat transfer : 
<p align=\"center\" style=\"font-style:italic;\">
g<sub>v</sub> = &beta;<sub>p</sub> ( p<sub>a</sub> - p<sub>s</sub> )
</p>
g<sub>v</sub> is the water flux density [kg/m<sup>2</sup>s], &beta;<sub>p</sub> is the water vapour transfer coefficient [kg/m<sup>2</sup>sPa],
p<sub>s</sub> is the water vapour pressure on the building component surface [Pa] and p<sub>a</sub> is the ambient water vapour pressure.
&beta;<sub>p</sub> is evaluated at 75E-9 for convection mass transfer in an outdoor surface. 
</li>
</ol>
<h4>References</h4>
<p>
Sparrow, E. M., J. W. Ramsey, and E. A. Mass. 1979. Effect of Finite Width on Heat Transfer 
and Fluid Flow about an Inclined Rectangular Plate. Journal of Heat Transfer, Vol. 101, p.
204.
</p>
<p>
Walton, G. N. 1981. Passive Solar Extension of the Building Loads Analysis and System
Thermodynamics (BLAST) Program, Technical Report, United States Army Construction
Engineering Research Laboratory, Champaign, IL.
</p>
<p>
Hartwig M. Kunzel. 1995. Simultaneous Heat and Moisture transport in Building Components. Franhofer Institute of Building Physics (ISBN 3-8167-4103-7).
</p>
</html>", revisions="<html>
<ul>

<li>
March 10 2014, by Antoine Hudzik:<br/>
First implementation.
</li>
</ul>
</html>"));
end ExteriorHM;
