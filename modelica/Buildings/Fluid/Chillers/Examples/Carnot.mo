within Buildings.Fluid.Chillers.Examples;
model Carnot "Test model for chiller based on Carnot efficiency"
  import Buildings;
 package Medium1 = Buildings.Media.ConstantPropertyLiquidWater
    "Medium model";
 package Medium2 = Buildings.Media.ConstantPropertyLiquidWater
    "Medium model";

  parameter Modelica.SIunits.Power P_nominal=10E3
    "Nominal compressor power (at y=1)";
  parameter Modelica.SIunits.TemperatureDifference dTEva_nominal=10
    "Temperature difference evaporator inlet-outlet";
  parameter Modelica.SIunits.TemperatureDifference dTCon_nominal=10
    "Temperature difference condenser outlet-inlet";
  parameter Real COPc = 3 "Chiller COP";
  parameter Modelica.SIunits.MassFlowRate m1_flow_nominal=
     P_nominal*(COPc+1)/dTEva_nominal/4200
    "Nominal mass flow rate at evaporator";
  parameter Modelica.SIunits.MassFlowRate m2_flow_nominal=
     m1_flow_nominal*COPc/(COPc+1) "Nominal mass flow rate at condenser";

  Buildings.Fluid.Chillers.Carnot chi(
    redeclare package Medium1 = Medium1,
    redeclare package Medium2 = Medium2,
    P_nominal=P_nominal,
    dTEva_nominal=dTEva_nominal,
    dTCon_nominal=dTCon_nominal,
    COP_nominal=COPc,
    use_eta_Carnot=true,
    etaCar=0.3,
    dp1_nominal=6000,
    dp2_nominal=6000,
    m1_flow_nominal=m1_flow_nominal,
    m2_flow_nominal=m2_flow_nominal) "Chiller model" 
    annotation (Placement(transformation(extent={{0,0},{20,20}})));
  Buildings.Fluid.Sources.MassFlowSource_T sou1(nPorts=1,
    redeclare package Medium = Medium1,
    use_T_in=true,
    m_flow=m1_flow_nominal,
    T=298.15) 
    annotation (Placement(transformation(extent={{-60,6},{-40,26}})));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics), Commands(file=
          "Carnot.mos" "run"));
  Buildings.Fluid.Sources.MassFlowSource_T sou2(nPorts=1,
    redeclare package Medium = Medium2,
    use_T_in=true,
    m_flow=m2_flow_nominal,
    T=291.15) 
    annotation (Placement(transformation(extent={{60,-6},{40,14}})));
  Buildings.Fluid.Sources.FixedBoundary sin1(nPorts=1, redeclare package Medium
      = Medium1)                                     annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={70,40})));
  Buildings.Fluid.Sources.FixedBoundary sin2(nPorts=1, redeclare package Medium
      = Medium2)                                     annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-50,-20})));
  Modelica.Blocks.Sources.Ramp uCom(
    height=-1,
    duration=60,
    offset=1,
    startTime=1800) "Compressor control signal" 
    annotation (Placement(transformation(extent={{-60,50},{-40,70}})));
  inner Modelica.Fluid.System system 
    annotation (Placement(transformation(extent={{-80,-80},{-60,-60}})));
  Modelica.Blocks.Sources.Ramp TEva_in(
    height=10,
    duration=60,
    offset=273.15 + 20,
    startTime=60) "Evaporator inlet temperature" 
    annotation (Placement(transformation(extent={{-90,10},{-70,30}})));
  Modelica.Blocks.Sources.Ramp TCon_in(
    height=10,
    duration=60,
    startTime=900,
    offset=273.15 + 10) "Condenser inlet temperature" 
    annotation (Placement(transformation(extent={{50,-40},{70,-20}})));
equation
  connect(sou1.ports[1], chi.port_a1)    annotation (Line(
      points={{-40,16},{0,16}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sou2.ports[1], chi.port_a2)    annotation (Line(
      points={{40,4},{20,4}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(chi.port_b1, sin1.ports[1])    annotation (Line(
      points={{20,16},{30,16},{30,40},{60,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sin2.ports[1], chi.port_b2)    annotation (Line(
      points={{-40,-20},{-10,-20},{-10,4},{0,4}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(TEva_in.y, sou1.T_in) annotation (Line(
      points={{-69,20},{-62,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TCon_in.y, sou2.T_in) annotation (Line(
      points={{71,-30},{80,-30},{80,8},{62,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(uCom.y, chi.y) annotation (Line(
      points={{-39,60},{-10,60},{-10,19},{-2,19}},
      color={0,0,127},
      smooth=Smooth.None));
end Carnot;