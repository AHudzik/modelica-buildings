within Buildings.HeatTransfer.Convection;
model InteriorHM "Model for an interior convective heat and mass transfer"
  extends BaseClasses.PartialConvectionHM;

    parameter Buildings.HeatTransfer.Types.InteriorConvection conMod=
    Buildings.HeatTransfer.Types.InteriorConvection.Fixed
    "Convective heat transfer model"
  annotation(Evaluate=true);

  parameter Boolean homotopyInitialization = true "= true, use homotopy method"
    annotation(Evaluate=true, Dialog(tab="Advanced"));

     parameter Real beta = 25E-9 "indoor water vapour coefficient (Kunzel)";
protected
  constant Modelica.SIunits.Temperature dT0 = 2
    "Initial temperature used in homotopy method";
   constant Real hv = 25000000 "latent heat of phase change";

equation
  if (conMod == Buildings.HeatTransfer.Types.InteriorConvection.Fixed) then
    q_flow = hFixed * dT;
    g_flow = hFixed * 7E-9 * dP + hv * m_flow;
  else

    // Even if hCon is a step function with a step at zero,
    // the product hCon*dT is differentiable at zero with
    // a continuous first derivative
    if homotopyInitialization then

      if isCeiling then
         q_flow = homotopy(actual=Buildings.HeatTransfer.Convection.Functions.HeatFlux.ceiling(dT=dT),
                    simplified=dT/dT0*Buildings.HeatTransfer.Convection.Functions.HeatFlux.ceiling(dT=dT0));
         g_flow = beta *dP + hv * m_flow;
      elseif isFloor then
         q_flow = homotopy(actual=Buildings.HeatTransfer.Convection.Functions.HeatFlux.floor(dT=dT),
                    simplified=dT/dT0*Buildings.HeatTransfer.Convection.Functions.HeatFlux.floor(dT=dT0));
         g_flow = beta*dP + hv * m_flow;
      else
         q_flow = homotopy(actual=Buildings.HeatTransfer.Convection.Functions.HeatFlux.wall(dT=dT),
                    simplified=dT/dT0*Buildings.HeatTransfer.Convection.Functions.HeatFlux.wall(dT=dT0));
         g_flow = beta*dP + hv * m_flow;
      end if;
    else
      if isCeiling then
         q_flow = Buildings.HeatTransfer.Convection.Functions.HeatFlux.ceiling(dT=dT);
         g_flow = beta*dP + hv * m_flow;
      elseif isFloor then
         q_flow = Buildings.HeatTransfer.Convection.Functions.HeatFlux.floor(dT=dT);
         g_flow = beta*dP + hv * m_flow;
      else
         q_flow = Buildings.HeatTransfer.Convection.Functions.HeatFlux.wall(dT=dT);
         g_flow = beta*dP + hv * m_flow;
      end if;

    end if;

   end if;

  annotation (
           Documentation(info="<html>
This is a model for a convective heat transfer for interior, room-facing surfaces.
The parameter <code>conMod</code> determines the model that is used to compute
the heat transfer coefficient:
<br/>

<ul>
<li><p>If <code>conMod=<a href=\"modelica://Buildings.HeatTransfer.Types.InteriorConvection\">
Buildings.HeatTransfer.Types.InteriorConvection.Fixed</a></code>, then
the convective heat transfer coefficient is set to the value specified by the parameter
<code>hFixed</code>.
</p>
</li>
<li>

If <code>conMod=<a href=\"modelica://Buildings.HeatTransfer.Types.InteriorConvection\">
Buildings.HeatTransfer.Types.InteriorConvection.Temperature</a></code>, then
the convective heat tranfer coefficient is a function of the temperature difference.
The convective heat flux is computed using
<br/>
<ul>
<li>
for floors the function 
<a href=\"modelica://Buildings.HeatTransfer.Convection.Functions.HeatFlux.floor\">
Buildings.HeatTransfer.Convection.Functions.HeatFlux.floor</a>
</li>
<li>
for ceilings the function
<a href=\"modelica://Buildings.HeatTransfer.Convection.Functions.HeatFlux.ceiling\">
Buildings.HeatTransfer.Convection.Functions.HeatFlux.ceiling</a>
</li>
<li>
for walls the function
<a href=\"modelica://Buildings.HeatTransfer.Convection.Functions.HeatFlux.wall\">
Buildings.HeatTransfer.Convection.Functions.HeatFlux.wall</a>
</li>
</li>
</ul>
</li>
</ul>
<p>The water vapour transfer is calculated in a manner similar to the heat transfer : 
<p align=\"center\" style=\"font-style:italic;\">
g<sub>v</sub> = &beta;<sub>p</sub> ( p<sub>a</sub> - p<sub>s</sub> )
</p>
g<sub>v</sub> is the water flux density [kg/m<sup>2</sup>s], &beta;<sub>p</sub> is the water vapour transfer coefficient [kg/m<sup>2</sup>sPa],
p<sub>s</sub> is the water vapour pressure on the building component surface [Pa] and p<sub>a</sub> is the ambient water vapour pressure.
&beta;<sub>p</sub> is evaluated at 25E-9 for convection mass transfer in an indoor surface (Kunzel). 
</p>
</html>", revisions="<html>
<ul>
<li>
April 2, 2011 by Michael Wetter:<br/>
Added <code>homotopy</code> operator.
</li>
<li>
March 10 2010, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
end InteriorHM;
