within Buildings.Rooms.Constructions.Examples;
model ExteriorWallHM "Test model for an exterior wall without a window"
  extends Modelica.Icons.Example;
protected
  BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam=
        "modelica://Buildings/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos",
    TDryBulSou=Buildings.BoundaryConditions.Types.DataSource.File,
    relHumSou=Buildings.BoundaryConditions.Types.DataSource.File)
    annotation (Placement(transformation(extent={{62,64},{86,92}})));
  BoundaryConditions.WeatherData.Bus weaBus
    annotation (Placement(transformation(extent={{112,30},{132,50}})));
  HeatTransfer.Interfaces.HeatMassPort_b opa_b1
    annotation (Placement(transformation(extent={{-24,-26},{-14,-16}})));
  HeatTransfer.Interfaces.HeatMassPort_a opa_a1
    annotation (Placement(transformation(extent={{-136,-28},{-126,-18}})));
  Utilities.Psychrometrics.X_pTphi x_pTphi
    annotation (Placement(transformation(extent={{86,-44},{66,-24}})));
  HeatTransfer.Sources.PrescribedHumidity prescribedHumidity
    annotation (Placement(transformation(extent={{42,-44},{22,-24}})));
  HeatTransfer.Sources.FixedHumidity fixedHumidity(X=0.006)
    annotation (Placement(transformation(extent={{-180,-60},{-160,-40}})));
  HeatTransfer.Sources.FixedTemperature fixedTemperature(T=294.15)
    annotation (Placement(transformation(extent={{-174,24},{-154,44}})));
  HeatTransfer.Sources.PrescribedTemperature prescribedTemperature
    annotation (Placement(transformation(extent={{40,20},{20,40}})));
  ConstructionHM conOpa(
    A=1,
    activatesuction=false,
    switch_w=2,
    switch_lamb=2,
    switch_dw=2,
    redeclare
      Buildings.HeatTransfer.Data.OpaqueConstructions.GypsumBoard50Concrete200HM
      layers,
    til(displayUnit="rad") = Buildings.HeatTransfer.Types.Tilt.Wall)
    annotation (Placement(transformation(extent={{-74,-54},{-28,-2}})));
equation

  connect(weaDat.weaBus, weaBus) annotation (Line(
      points={{86,78},{122,78},{122,40}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(weaBus.pAtm, x_pTphi.p_in) annotation (Line(
      points={{122,40},{106,40},{106,-28},{88,-28}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(weaBus.relHum, x_pTphi.phi) annotation (Line(
      points={{122,40},{106,40},{106,-40},{88,-40}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(x_pTphi.X[1], prescribedHumidity.X) annotation (Line(
      points={{65,-34},{44,-34}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(prescribedHumidity.massPort, opa_b1.massPort) annotation (Line(
      points={{22,-33.8},{4,-33.8},{4,-21},{-19,-21}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedTemperature.port, opa_a1.heatPort) annotation (Line(
      points={{-154,34},{-146,34},{-146,-23},{-131,-23}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedHumidity.massPort, opa_a1.massPort) annotation (Line(
      points={{-160,-50},{-146,-50},{-146,-23},{-131,-23}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(weaBus.TDryBul, prescribedTemperature.T) annotation (Line(
      points={{122,40},{82,40},{82,30},{42,30}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(weaBus.TDryBul, x_pTphi.T) annotation (Line(
      points={{122,40},{106,40},{106,-34},{88,-34}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(prescribedTemperature.port, opa_b1.heatPort) annotation (Line(
      points={{20,30},{2,30},{2,-21},{-19,-21}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(opa_b1, conOpa.opa_b) annotation (Line(
      points={{-19,-21},{-24,-21},{-24,-10.58},{-27.8467,-10.58}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(opa_a1, conOpa.opa_a) annotation (Line(
      points={{-131,-23},{-130,-23},{-130,-10.4067},{-74.23,-10.4067}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(opa_b1, opa_b1) annotation (Line(
      points={{-19,-21},{-22.5,-21},{-22.5,-21},{-19,-21}},
      color={127,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}},
          preserveAspectRatio=false), graphics), Icon(coordinateSystem(extent={{-100,
            -100},{100,100}})));
end ExteriorWallHM;
