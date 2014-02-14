within Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses;
function singleUTubeResistances
  "Thermal resistances for single U-tube, according to Bauer et al. (2011)"

  // Geometry of the borehole
  input Modelica.SIunits.Height hSeg "Height of the element";
  input Modelica.SIunits.Radius rBor "Radius of the borehole";
  // Geometry of the pipe
  input Modelica.SIunits.Radius rTub "Radius of the tube";
  input Modelica.SIunits.Length eTub "Thickness of the tubes";
  input Modelica.SIunits.Length sha
    "Shank spacing, defined as the distance between the center of a pipe and the center of the borehole";

  // Thermal properties
  input Modelica.SIunits.ThermalConductivity kFil
    "Thermal conductivity of the grout";
  input Modelica.SIunits.ThermalConductivity kSoi
    "Thermal conductivity of the soi";
  input Modelica.SIunits.ThermalConductivity kTub
    "Thermal conductivity of the tube";

  // Outputs
  output Modelica.SIunits.ThermalResistance Rgb
    "Thermal resistance between grout zone and borehole wall";
  output Modelica.SIunits.ThermalResistance Rgg
    "Thermal resistance between the two grout zones";
  output Modelica.SIunits.ThermalResistance RCondGro
    "Thermal resistance between: pipe wall to capacity in grout";
  output Real x "Capacity location";

protected
  Boolean test=false "thermodynamic test for R and x value";

  Modelica.SIunits.ThermalResistance Rg
    "Thermal resistance between outer borehole wall and one tube";
  Modelica.SIunits.ThermalResistance Rar
    "Thermal resistance between the two pipe outer walls";
  Modelica.SIunits.ThermalResistance RCondPipe
    "Thermal resistance of the pipe wall";

  Real Rb
    "Fluid-to-grout resistance, as defined by Hellstroem. Resistance from the fluid in the pipe to the borehole wall";
  Real Ra
    "Grout-to-grout resistance (2D) as defined by Hellstroem. Interaction between the different grout part";

  // Help variables
  Real sigma "Help variable as defined by Hellstroem";
  Real beta "Help variable as defined by Hellstroem";
  Real R_1delta_LS
    "One leg of the triangle resistance network, corresponding to the line source solution";
  Real R_1delta_MP
    "One leg of the triangle resistance network, corresponding to the multipole solution";
  Real Ra_LS
    "Grout-to-grout resistance calculated with the line-source approximation";

  Integer i=1 "Loop counter";

algorithm
  // ********** Rb and Ra from multipole **********
  // Help variables
  RCondPipe :=Modelica.Math.log((rTub + eTub)/rTub)/(2*Modelica.Constants.pi*hSeg*kTub);
  sigma :=(kFil - kSoi)/(kFil + kSoi);
  R_1delta_LS :=1/(2*Modelica.Constants.pi*kFil)*(log(rBor/(rTub + eTub)) + log(rBor/(2*sha)) +
    sigma*log(rBor^4/(rBor^4 - sha^4)));
  R_1delta_MP :=R_1delta_LS - 1/(2*Modelica.Constants.pi*kFil)*((rTub + eTub)^2/
    (4*sha^2)*(1 - sigma*4*sha^4/(rBor^4 - sha^4))^2)/((1 + beta)/(1 - beta) + (
    rTub + eTub)^2/(4*sha^2)*(1 + sigma*16*sha^4*rBor^4/(rBor^4 - sha^4)^2));
  Ra_LS      :=1/(Modelica.Constants.pi*kFil)*(log(2*sha/rTub) + sigma*log((
    rBor^2 + sha^2)/(rBor^2 - sha^2)));

  //Rb and Ra
  beta :=2*Modelica.Constants.pi*kFil*RCondPipe;
  Rb :=R_1delta_MP/2;
  Ra :=Ra_LS - 1/(Modelica.Constants.pi*kFil)*(rTub^2/(4*sha^2)*(1 + sigma*
    4*rBor^4*sha^2/(rBor^4 - sha^4))/((1 + beta)/(1 - beta) - rTub^2/(4*sha^2) +
    sigma*2*rTub^2*rBor^2*(rBor^4 + sha^4)/(rBor^4 - sha^4)^2));

  //Conversion of Rb (resp. Ra) to Rg (resp. Rar) of Bauer:
  Rg  :=2*Rb/hSeg;
  Rar :=Ra/hSeg;

/* **************** Simplification of Bauer for single U-tube ************************
  //Thermal resistance between: Outer wall and one tube
     Rg := Modelica.Math.acosh((rBor^2 + (rTub + eTub)^2 - sha^2)/(2*rBor*(rTub +
       eTub)))/(2*Modelica.Constants.pi*hSeg*kFil)*(1.601 - 0.888*sha/rBor);

  //Thermal resistance between: The two pipe outer walls
  Rar := Modelica.Math.acosh((2*sha^2 - (rTub + eTub)^2)/(rTub + eTub)^2)/(2*
       Modelica.Constants.pi*hSeg*kFil);
*************************************************************************************** */

  // ********** Resistances and capacity location according to Bauer **********
  while test == false and i <= 10 loop
    // Capacity location (with correction factor in case that the test is negative)
    x := Modelica.Math.log(sqrt(rBor^2 + 2*(rTub + eTub)^2)/(2*(rTub + eTub)))/
      Modelica.Math.log(rBor/(sqrt(2)*(rTub + eTub)))*((15 - i + 1)/15);

    //Thermal resistance between the grout zone and bore hole wall
    Rgb := (1 - x)*Rg;

    //Thermal resistance between the two grout zones
    Rgg := 2*Rgb*(Rar - 2*x*Rg)/(2*Rgb - Rar + 2*x*Rg);

    // Thermodynamic test to check if negative R values make sense. If not, decrease x-value.
    // fixme: the implemented is only for single U-tube BHE's.
    test := ((1/Rgg + 1/2/Rgb)^(-1) > 0);
    i := i + 1;
  end while;
  //Conduction resistance in grout from pipe wall to capacity in grout
  RCondGro := x*Rg + RCondPipe;

  annotation (Diagram(graphics), Documentation(info="<html>
<p>
This model computes the different thermal resistances present in a 
single-U-tube borehole using the method of Bauer et al. [1]. 
It also computes the fluid-to-ground thermal resistance <i>Rb</i> 
and the grout-to-grout thermal resistance <i>Ra</i> as defined by 
Hellstroem [2] using the multipole method.
</p>
<p>
The figure below shows the thermal network set up by Bauer et al.
</p>
<p align=\"center\">
<img src=\"modelica://Buildings 1.5/Buildings/Resources/Images/Fluid/HeatExchangers/Boreholes/BaseClassesBauer_singleUTube_small.png\" alt=\"image\"/> <!-- fixme -->
</p>
<p>
The different resistances are calculated with following equations:
</p>
<p>
Grout zone and bore hole wall:<i> R<sub>gb</sub><sup>1U</sup> =  (  1 - x<sup>1U</sup>  )  R<sub>g</sub><sup>1U</sup></i> </p>
<p>
Thermal resistance between the two grout zones: 
<i>R<sub>gg</sub><sup>1U</sup> := 2 R<sub>gb</sub><sup>1U</sup>  ( R<sub>ar</sub><sup>1U</sup> - 2 x<sup>1U</sup> R<sub>g</sub><sup>1U</sup> ) / ( 2 R<sub>gb</sub><sup>1U</sup> - R<sub>ar</sub><sup>1U</sup> + 2 x<sup>1U</sup> R<sub>g</sub><sup>1U</sup> ) 
</i>
</p>
<p>
Thermal resistance between the pipe wall to capacity in grout:
<i> 
RCondGro = x<sup>1U</sup> R<sub>g</sub><sup>1U</sup> + log (  ( rTub + eTub ) /rTub ) / ( 2 &pi; hSeg kTub )  
</i>
</p>
<p>
Capacity location:  <-- fixme -->
<i> 
x<sup>1U</sup> =log ( &radic;<span style=\"text-decoration:overline;\">rBor<sup>2</sup> 
+ 2  ( rTub + eTub ) <sup>2</sup></span>/ ( 2  ( rTub + eTub )  )  ) / 
log ( rBor/ ( &radic;<span style=\"text-decoration:overline;\">2</span>  ( rTub + eTub )  )  )  
</i>
</p>
<p>
Thermal resistance between outer borehole wall and one tube: 
<i>
R<sub>g</sub><sup>1U</sup> =2 Rb/hSeg </i>
</p>
<p>
Thermal resistance between the two pipe outer walls:
<i>
R<sub>ar</sub><sup>1U</sup> :=Ra/hSeg
</i> 
</p>
<p>
The fluid-to-ground thermal resistance 
<i>Rb</i> 
and the grout-to-grout resistance <i>Ra</i> are calculated with the multipole method 
(Hellstroem (1991)) shown below.
</p>
<p>
<i>
R_b =1/ ( 4 &pi; kFil )   ( log ( rBor/ ( rTub + eTub )  )  + log ( rBor/ ( 2 sha )  )  +
    &sigma; log ( rBor<sup>4</sup>/ ( rBor<sup>4</sup> - sha<sup>4</sup> )  )  )  - 1/ ( 4 &pi; kFil )   (  ( rTub + eTub ) <sup>2</sup>/
     ( 4 sha<sup>2</sup> )   ( 1 - &sigma; 4 sha<sup>4</sup>/ ( rBor<sup>4</sup> - sha<sup>4</sup> )  ) <sup>2</sup> ) / (  ( 1 + &beta; ) / ( 1 - &beta; )  +  ( 
    rTub + eTub ) <sup>2</sup>/ ( 4 sha<sup>2</sup> )   ( 1 + &sigma; 16 sha<sup>4</sup> rBor<sup>4</sup>/ ( rBor<sup>4</sup> - sha<sup>4</sup> ) <sup>2</sup> )  ) 
</i>
</p>
<p>
<i>
R_a = 1/ ( &pi; kFil )   ( log ( 2 sha/rTub )  + &sigma; log ((
    rBor<sup>2</sup> + sha<sup>2</sup> ) / ( rBor<sup>2</sup> - sha<sup>2</sup> )  )  )  - 1/ ( &pi; kFil )   ( rTub<sup>2</sup>/ ( 4 sha<sup>2</sup> )   ( 1 + &sigma; 
    4 rBor<sup>4</sup> sha<sup>2</sup>/ ( rBor<sup>4</sup> - sha<sup>4</sup> )  ) / (  ( 1 + &beta; ) / ( 1 - &beta; )  - rTub<sup>2</sup>/ ( 4 sha<sup>2</sup> )  +
    &sigma; 2 rTub<sup>2</sup> rBor<sup>2</sup>  ( rBor<sup>4</sup> + sha<sup>4</sup> ) / ( rBor<sup>4</sup> - sha<sup>4</sup> ) <sup>2</sup> )  ) 
</i>
</p>
<p> 
with 
<i>
&sigma; = ( kFil - kSoi ) / ( kFil + kSoi ) </i> and  <i> &beta; :=2 &pi; kFil RCondPipe
</i>
</p>
<p>
where
<i>kFil</i> and <i>kSoi</i> are the conductivity of the filling material 
and of the ground respectively,
<i>rTub+eTub</i> and <i>rBor</i> are the pipe and the borehole outside radius and 
<i>sha</i> is the shank spacing, which is equal to distance between center of borehole to center of pipe).
</p>
<p>
<h4>References</h4>
</p>
<p>
G. Hellstr&ouml;m. <i>Ground heat storage: thermal analyses of duct storage systems (Theory)</i>. 
Dept. of Mathematical Physics, University of Lund, Sweden, 1991.
</p>
<p>D. Bauer, W. Heidemann, H. M&uuml;ller-Steinhagen, and H.-J. G. Diersch. 
<i>Thermal resistance and capacity models for borehole heat exchangers</i>. 
International Journal Of Energy Research, 35:312&ndash;320, 2010.
</p>
</html>", revisions="<html>
<p>
<ul>
<li>
February 13, 2014 by Damien Picard:<br/>
Edit documentation and add formule for beta.
</li>
<li>
February 12, 2014, by Damien Picard:<br/>
Remove the flow dependency of the resistances, as this function calculates the conduction resistances only.
</li>
<li>
January 24, 2014, by Michael Wetter:<br/>
Revised implementation.</li>
<li>
January 23, 2014, by Damien Picard:<br/>
First implementation.
</li>
</ul>
</p>
</html>"));
end singleUTubeResistances;