within Buildings.Rooms.Constructions.Examples;
model ExteriorWallHM "Test model for an exterior wall without a window"
  import Buildings;
  extends Modelica.Icons.Example;

  ConstructionHM conOpa(
    activatesuction=false,
    redeclare Buildings.HeatTransfer.Data.OpaqueConstructions.Concrete200HM
      layers,
    til=Buildings.HeatTransfer.Types.Tilt.Wall,
    A=3*10,
    switch_w=1,
    switch_lamb=1,
    switch_dw=1)
    annotation (Placement(transformation(extent={{-66,-62},{4,12}})));

protected
  BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam=
        "modelica://Buildings/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos",
    TDryBulSou=Buildings.BoundaryConditions.Types.DataSource.File,
    relHumSou=Buildings.BoundaryConditions.Types.DataSource.File)
    annotation (Placement(transformation(extent={{92,-80},{116,-52}})));
  HeatTransfer.Sources.FixedHumidity fixedHumidity(X=0.006)
    annotation (Placement(transformation(extent={{-98,-40},{-78,-20}})));
  HeatTransfer.Sources.FixedTemperature fixedTemperature(T=294.15)
    annotation (Placement(transformation(extent={{-100,20},{-80,40}})));

public
  Buildings.Rooms.BaseClasses.ExteriorBoundaryConditionsHM
    exteriorBoundaryConditionsHM(
    conPar=conPar,
    nCon=1,
    lat(displayUnit="rad") = 0.73268921998722,
    linearizeRadiation=false)
    annotation (Placement(transformation(extent={{14,-50},{78,12}})));
  Buildings.HeatTransfer.Data.OpaqueConstructions.Concrete200HM concrete200HM
    annotation (Placement(transformation(extent={{-40,68},{-20,88}})));

  parameter Buildings.Rooms.BaseClasses.ParameterConstruction conPar[1](
    each til=Buildings.HeatTransfer.Types.Tilt.Wall,
    each azi=0,
    each A=3*10,
    redeclare Buildings.HeatTransfer.Data.OpaqueConstructions.Concrete200
      layers) "Data for construction"
    annotation (Placement(transformation(extent={{-10,68},{10,88}})));
  Buildings.HeatTransfer.Interfaces.HeatMassPort_a opa_a1
    annotation (Placement(transformation(extent={{-88,-10},{-68,10}})));
equation

  connect(conOpa.opa_b, exteriorBoundaryConditionsHM.fluid1) annotation (Line(
      points={{4.23333,-0.21},{8.11333,-0.21},{8.11333,-0.4},{14,-0.4}},
      color={127,0,127},
      smooth=Smooth.None));
  connect(weaDat.weaBus, exteriorBoundaryConditionsHM.weaBus) annotation (Line(
      points={{116,-66},{78,-66},{78,-16.83},{69.7867,-16.83}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(conOpa.opa_a, opa_a1) annotation (Line(
      points={{-66.35,0.0366667},{-72.175,0.0366667},{-72.175,0},{-78,0}},
      color={0,0,0},
      pattern=LinePattern.None,
      smooth=Smooth.None));
  connect(fixedHumidity.massPort, opa_a1.massPort) annotation (Line(
      points={{-78,-30},{-78,0},{-78,0}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(fixedTemperature.port, opa_a1.heatPort) annotation (Line(
      points={{-80,30},{-80,0},{-78,0}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}},
          preserveAspectRatio=false), graphics), Icon(coordinateSystem(extent={{-100,
            -100},{100,100}})));
end ExteriorWallHM;
