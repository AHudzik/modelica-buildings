within Buildings.Rooms.Examples.TestConditionalConstructionsFFD;
model ShoeBox "A shoebox room model with only walls"
  extends Modelica.Icons.Example;
  package MediumA = Buildings.Media.GasesConstantDensity.MoistAirUnsaturated
    "Medium model";

  inner Modelica.Fluid.System system
    annotation (Placement(transformation(extent={{-72,-32},{-52,-12}})));

  parameter
    Buildings.HeatTransfer.Data.OpaqueConstructions.Insulation100Concrete200
    matLayExt "Construction material for exterior walls"
    annotation (Placement(transformation(extent={{-60,140},{-40,160}})));

  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Brick120 matLayPar
    "Construction material for partition walls"
    annotation (Placement(transformation(extent={{-20,140},{0,160}})));

  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic matLayRoo(
      material={HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
        HeatTransfer.Data.Solids.Concrete(x=0.2)}, final nLay=2)
    "Construction material for roof"
    annotation (Placement(transformation(extent={{20,140},{40,160}})));

  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic matLayFlo(
      material={HeatTransfer.Data.Solids.Concrete(x=0.2),
        HeatTransfer.Data.Solids.InsulationBoard(x=0.15),
        HeatTransfer.Data.Solids.Concrete(x=0.05)}, final nLay=3)
    "Construction material for floor"
    annotation (Placement(transformation(extent={{60,140},{80,160}})));

  parameter Buildings.HeatTransfer.Data.GlazingSystems.DoubleClearAir13Clear
    glaSys(
    UFra=2,
    shade=Buildings.HeatTransfer.Data.Shades.Gray(),
    haveInteriorShade=false,
    haveExteriorShade=false) "Data record for the glazing system"
    annotation (Placement(transformation(extent={{100,140},{120,160}})));

  parameter Integer nConExtWin=1 "Number of constructions with a window";
  parameter Integer nConBou=1
    "Number of surface that are connected to constructions that are modeled inside the room";
  parameter Integer nSurBou=1
    "Number of surface that are connected to the room air volume";
  parameter String cfdFilNam "CFD input file name";
  Buildings.Rooms.FFD roo(
    redeclare package Medium = MediumA,
    nConBou=nConBou,
    nSurBou=nSurBou,
    sensorName={"Occupied zone air temperature","Velocity"},
    linearizeRadiation=false,
    useFFD=true,
    startTime=0,
    samplePeriod=1,
    nConPar=0,
    nConExt=4,
    nConExtWin=0,
    AFlo=1*1,
    hRoo=1,
    datConBou(
      name={"Floor"},
      layers={matLayFlo},
      each A=1*1,
      each til=Buildings.HeatTransfer.Types.Tilt.Floor),
    datConExt(
      name={"Ceiling","West Wall","North Wall","South Wall"},
      layers={matLayRoo,matLayExt,matLayExt,matLayExt},
      A={1*1,1*1,1*1,1*1},
      til={Buildings.HeatTransfer.Types.Tilt.Ceiling,Buildings.HeatTransfer.Types.Tilt.Wall,
          Buildings.HeatTransfer.Types.Tilt.Wall,Buildings.HeatTransfer.Types.Tilt.Wall},
      azi={Buildings.HeatTransfer.Types.Azimuth.S,Buildings.HeatTransfer.Types.Azimuth.W,
          Buildings.HeatTransfer.Types.Azimuth.N,Buildings.HeatTransfer.Types.Azimuth.S}),
    surBou(
      name={"East Wall"},
      each A=1*1,
      each absIR=0.9,
      each absSol=0.9,
      each til=Buildings.HeatTransfer.Types.Tilt.Wall),
    cfdFilNam = cfdFilNam,
    lat=0.73268921998722,
    cfddFilNam="Buildings/Resources/Data/Rooms/FFD/ShoeBox.ffd") "Room model"
    annotation (Placement(transformation(extent={{46,20},{86,60}})));

  Modelica.Blocks.Sources.Constant qConGai_flow(k=0) "Convective heat gain"
    annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
  Modelica.Blocks.Sources.Constant qRadGai_flow(k=0) "Radiative heat gain"
    annotation (Placement(transformation(extent={{-60,80},{-40,100}})));
  Modelica.Blocks.Routing.Multiplex3 multiplex3_1
    annotation (Placement(transformation(extent={{-20,40},{0,60}})));
  Modelica.Blocks.Sources.Constant qLatGai_flow(k=0) "Latent heat gain"
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
    filNam="Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos",
    TDryBulSou=Buildings.BoundaryConditions.Types.DataSource.Parameter,
    TDryBul=293.15)
    annotation (Placement(transformation(extent={{160,140},{180,160}})));

  Modelica.Blocks.Sources.Constant uSha(k=0)
    "Control signal for the shading device"
    annotation (Placement(transformation(extent={{-20,90},{0,110}})));
  Modelica.Blocks.Routing.Replicator replicator(nout=max(1, nConExtWin))
    annotation (Placement(transformation(extent={{10,90},{30,110}})));
  Buildings.HeatTransfer.Sources.FixedTemperature TSoi[nConBou](each T=283.15)
    "Boundary condition for construction" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={110,-10})));
  Buildings.HeatTransfer.Sources.FixedTemperature TBou[nSurBou](each T=288.15)
    "Boundary condition for construction" annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={150,-50})));
  HeatTransfer.Conduction.MultiLayer conOut[nSurBou](redeclare
      Buildings.HeatTransfer.Data.OpaqueConstructions.Brick120 layers, each A=6
        *4) "Construction that is modeled outside of room"
    annotation (Placement(transformation(extent={{100,-60},{120,-40}})));

equation
  connect(qRadGai_flow.y, multiplex3_1.u1[1]) annotation (Line(
      points={{-39,90},{-32,90},{-32,57},{-22,57}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(qConGai_flow.y, multiplex3_1.u2[1]) annotation (Line(
      points={{-39,50},{-22,50}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(qLatGai_flow.y, multiplex3_1.u3[1]) annotation (Line(
      points={{-39,10},{-32,10},{-32,43},{-22,43}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex3_1.y, roo.qGai_flow) annotation (Line(
      points={{1,50},{20,50},{20,48},{44,48}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(weaDat.weaBus, roo.weaBus) annotation (Line(
      points={{180,150},{190,150},{190,57.9},{83.9,57.9}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(uSha.y, replicator.u) annotation (Line(
      points={{1,100},{8,100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TSoi.port, roo.surf_conBou) annotation (Line(
      points={{100,-10},{72,-10},{72,24}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(TBou.port, conOut.port_b) annotation (Line(
      points={{140,-50},{120,-50}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(roo.surf_surBou, conOut.port_a) annotation (Line(
      points={{62.2,26},{62,26},{62,-50},{100,-50}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(roo.uSha, replicator.y) annotation (Line(
      points={{44,56},{40,56},{40,100},{31,100}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            200,200}}), graphics),
    __Dymola_Commands(file=
          "modelica://Buildings/Resources/Scripts/Dymola/Rooms/Examples/TestConditionalConstructionsFFD/FreeResponse.mos"
        "Simulate and plot"),
    Documentation(info="<html>
<p>
This model tests whether 
<a href=\"modelica://Buildings.Rooms.FFD\">
Buildings.Rooms.FFD</a>
can conduct cosimulation with FFD program with only the walls. 
The dimension of the room is 1m x 1m x 1m.
</p>
<p>
<b>Note:</b> This model has an unrealistic geometry as it is used only 
to test the correctness of cosimulaton with walls.
</p>
</html>", revisions="<html>
<ul>
<li>
August 13, 2013, by Wangda Zuo:<br/>
First implementation.
</li>
</ul>
</html>"),
    experiment(
      StopTime=3,
      __Dymola_fixedstepsize=0.1,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput);
end ShoeBox;
