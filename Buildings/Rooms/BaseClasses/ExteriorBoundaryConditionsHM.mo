within Buildings.Rooms.BaseClasses;
model ExteriorBoundaryConditionsHM
  "Model for convection and radiation bounary condition of exterior constructions"
  parameter Integer nCon(min=1) "Number of exterior constructions"
  annotation (Dialog(group="Exterior constructions"));
  parameter Modelica.SIunits.Angle lat "Latitude";

  parameter Boolean linearizeRadiation
    "Set to true to linearize emissive power";

  replaceable parameter ParameterConstruction conPar[nCon] constrainedby
    ParameterConstruction "Records for construction"
    annotation (Placement(transformation(extent={{174,-214},{194,-194}})));

  parameter Buildings.HeatTransfer.Types.ExteriorConvection conMod=
  Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind
    "Convective heat transfer model for opaque part of the constructions"
    annotation (Dialog(group="Convective heat transfer"));
  parameter Modelica.SIunits.CoefficientOfHeatTransfer hFixed=10.0
    "Constant convection coefficient for opaque part of the constructions"
    annotation (Dialog(group="Convective heat transfer",
                       enable=(conMod == Buildings.HeatTransfer.Types.ExteriorConvection.Fixed)));

  // The convection coefficients are not final to allow a user to individually
  // assign them.
  // We reassign the tilt since a roof has been declared in the room model as the
  // ceiling (of the room)

  SkyRadiationExchange skyRadExc(
    final n=nCon,
    final A=AOpa,
    final absIR=conPar[:].layers.absIR_a,
    vieFacSky={(Modelica.Constants.pi - conPar[i].til)./Modelica.Constants.pi for i in 1:nCon})
    "Infrared radiative heat exchange with sky"
    annotation (Placement(transformation(extent={{-140,240},{-180,280}})));
  BoundaryConditions.WeatherData.Bus weaBus
    annotation (Placement(transformation(extent={{234,32},{254,52}}),
        iconTransformation(extent={{192,-10},{254,52}})));

  BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil[
            nCon](
    each final lat=lat,
    final til=conPar[:].til,
    final azi=conPar[:].azi) "Direct solar irradiation on the surface"
    annotation (Placement(transformation(extent={{220,120},{200,140}})));
  BoundaryConditions.SolarIrradiation.DiffusePerez HDifTil[nCon](
    each final lat=lat,
    final til=conPar[:].til,
    final azi=conPar[:].azi) "Diffuse solar irradiation"
    annotation (Placement(transformation(extent={{220,80},{200,100}})));
  Modelica.Blocks.Math.Add HTotConExt[nCon](
    final k1=conPar[:].layers.absSol_a .* AOpa,
    final k2=conPar[:].layers.absSol_a .* AOpa) "Total solar irradiation"
    annotation (Placement(transformation(extent={{40,100},{20,120}})));
  Buildings.HeatTransfer.Sources.PrescribedHeatFlow solHeaGaiConExt[nCon]
    "Total solar heat gain of the surface"
    annotation (Placement(transformation(extent={{0,100},{-20,120}})));

protected
  parameter Modelica.SIunits.Area AOpa[nCon]=conPar[:].A
    "Area of opaque construction";

  Buildings.HeatTransfer.Sources.PrescribedTemperature TAirConExt[
    nCon] "Outside air temperature for exterior constructions"
    annotation (Placement(transformation(extent={{8,160},{-32,200}})));
  Modelica.Blocks.Routing.Replicator repConExt(nout=nCon) "Signal replicator"
    annotation (Placement(transformation(extent={{100,170},{80,190}})));

  Modelica.Blocks.Routing.Replicator repConExt1(
                                               nout=nCon) "Signal replicator"
    annotation (Placement(transformation(extent={{130,200},{110,220}})));
  Modelica.Blocks.Routing.Replicator repConExt2(
                                               nout=nCon) "Signal replicator"
    annotation (Placement(transformation(extent={{180,220},{160,240}})));

public
  HeatTransfer.Convection.ExteriorHM exteriorHM(
    til=conPar[1].til,
    azi=conPar[1].azi,
    A=AOpa[1])
    annotation (Placement(transformation(extent={{-120,160},{-160,200}})));
  Utilities.Psychrometrics.X_pTphi x_pTphi
    annotation (Placement(transformation(extent={{72,-2},{30,40}})));
  HeatTransfer.Interfaces.HeatMassPort_a solid1
    annotation (Placement(transformation(extent={{-152,-246},{-132,-226}}),
        iconTransformation(extent={{-152,-246},{-132,-226}})));
  HeatTransfer.Interfaces.HeatMassPort_b fluid1
    annotation (Placement(transformation(extent={{-310,170},{-290,190}})));
  HeatTransfer.Sources.PrescribedHumidity prescribedHumidity
    annotation (Placement(transformation(extent={{0,0},{-40,40}})));
equation

  connect(repConExt.y, TAirConExt.T) annotation (Line(
      points={{79,180},{12,180}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(repConExt.u, weaBus.TDryBul) annotation (Line(
      points={{102,180},{244,180},{244,42}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(skyRadExc.TOut, weaBus.TDryBul) annotation (Line(
      points={{-136,252},{244,252},{244,42}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(skyRadExc.TBlaSky, weaBus.TBlaSky) annotation (Line(
      points={{-136,268},{244,268},{244,42}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  for i in 1:nCon loop
  connect(weaBus, HDirTil[i].weaBus) annotation (Line(
      points={{244,42},{244,130},{220,130}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(HDifTil[i].weaBus, weaBus) annotation (Line(
      points={{220,90},{244,90},{244,42}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
   end for;
  connect(HTotConExt.y, solHeaGaiConExt.Q_flow) annotation (Line(
      points={{19,110},{5.55112e-16,110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(HDirTil.H, HTotConExt.u1) annotation (Line(
      points={{199,130},{60,130},{60,116},{42,116}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(HDifTil.H, HTotConExt.u2) annotation (Line(
      points={{199,90},{60,90},{60,104},{42,104}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(repConExt2.u, weaBus.winDir) annotation (Line(
      points={{182,230},{244,230},{244,42}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(repConExt1.u, weaBus.winSpe) annotation (Line(
      points={{132,210},{244,210},{244,42}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(weaBus.pAtm, x_pTphi.p_in) annotation (Line(
      points={{244,42},{184,42},{184,31.6},{76.2,31.6}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(weaBus.TDryBul, x_pTphi.T) annotation (Line(
      points={{244,42},{184,42},{184,19},{76.2,19}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(weaBus.relHum, x_pTphi.phi) annotation (Line(
      points={{244,42},{184,42},{184,6.4},{76.2,6.4}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(exteriorHM.solid, solid1) annotation (Line(
      points={{-120,181.2},{-98,181.2},{-98,-236},{-142,-236}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(exteriorHM.fluid, fluid1) annotation (Line(
      points={{-160,180},{-300,180}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(TAirConExt[1].port, solid1.heatPort) annotation (Line(
      points={{-32,180},{-88,180},{-88,-236},{-142,-236}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(repConExt1.y[1], exteriorHM.v) annotation (Line(
      points={{109,210},{-80,210},{-80,200},{-116,200}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(repConExt2.y[1], exteriorHM.dir) annotation (Line(
      points={{159,230},{-92,230},{-92,190},{-116,190}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(x_pTphi.X[1], prescribedHumidity.X) annotation (Line(
      points={{27.9,19},{15.95,19},{15.95,20},{4,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(prescribedHumidity.massPort, solid1.massPort) annotation (Line(
      points={{-40,20.4},{-142,20.4},{-142,-236}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(solHeaGaiConExt[1].port, fluid1.heatPort) annotation (Line(
      points={{-20,110},{-300,110},{-300,180}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(skyRadExc.port[1], fluid1.heatPort) annotation (Line(
      points={{-180,260},{-238,260},{-238,180},{-300,180}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-300,
            -300},{300,300}},
        initialScale=0.1), graphics),
                          Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-300,-300},{300,300}},
        initialScale=0.1), graphics={
        Rectangle(
          extent={{-160,280},{280,-250}},
          fillColor={230,243,255},
          fillPattern=FillPattern.Solid,
          lineColor={0,0,0}),  Ellipse(
          extent={{164,262},{270,162}},
          lineColor={255,255,0},
          fillColor={255,213,170},
          fillPattern=FillPattern.Sphere),
        Rectangle(
          extent={{-220,280},{-160,-280}},
          lineColor={0,0,0},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-160,-250},{280,-280}},
          lineColor={0,0,0},
          fillColor={0,127,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-168,346},{212,280}},
          lineColor={0,0,255},
          textString="%name"),
        Line(
          points={{-136,268},{-136,-234}},
          color={0,0,255},
          pattern=LinePattern.Dot,
          smooth=Smooth.None),
        Line(
          points={{-148,268},{-148,-230}},
          color={0,0,127},
          pattern=LinePattern.Dot,
          smooth=Smooth.None)}),
        Documentation(info="<html>
This model computes the boundary conditions for the outside-facing surface of
opaque constructions.
<p>
The model computes the infrared, solar, and convective heat and mass exchange
between these surfaces and the exterior temperature and the sky temperature.
Input into this model are weather data that may be obtained from
<a href=\"modelica://Buildings.BoundaryConditions.WeatherData\">
Buildings.BoundaryConditions.WeatherData</a>.
</p>
<p>
In this model, the solar radiation data are converted from horizontal irradiation to
irradiation on tilted surfaces using models from the package
<a href=\"modelica://Buildings.BoundaryConditions.SolarIrradiation\">
Buildings.BoundaryConditions.SolarIrradiation</a>.
The convective heat transfer between the exterior surface of the opaque constructions
is computed using
<a href=\"modelica://Buildings.HeatTransfer.Convection\">
Buildings.HeatTransfer.Convection</a>.
</p>
<p>
The heat transfer of windows are not computed in this model. They are implemented in
<a href=\"modelica://Buildings.Rooms.BaseClasses.ExteriorBoundaryConditionsWithWindow\">
Buildings.Rooms.BaseClasses.ExteriorBoundaryConditionsWithWindow</a>.
</p>
</html>",
        revisions="<html>
<ul>

<li>
July 28, 2011, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>
</html>"));
end ExteriorBoundaryConditionsHM;
