within Buildings.Fluid.Sources.Examples;
model TraceSubstancesFlowSource
  extends Modelica.Icons.Example;
  package Medium = Buildings.Media.GasesPTDecoupled.SimpleAir(extraPropertiesNames={"CO2"});

  MixingVolumes.MixingVolume vol(
    redeclare package Medium = Medium,
    V=100,
    m_flow_nominal=1,
    nPorts=2,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial) "Mixing volume"
                          annotation (Placement(transformation(extent={{100,108},
            {120,128}},rotation=0)));
  Sources.TraceSubstancesFlowSource sou(redeclare package Medium = Medium,
      use_m_flow_in=true,
    nPorts=1)
    annotation (Placement(transformation(extent={{-46,98},{-26,118}},rotation=0)));
  Modelica.Blocks.Sources.Step step(          startTime=0.5,
    height=-2,
    offset=2)
    annotation (Placement(transformation(extent={{-92,30},{-72,50}},  rotation=
            0)));
  FixedResistances.FixedResistanceDpM res(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=1)
    "Resistance, used to check if species are transported between ports"
    annotation (Placement(transformation(extent={{60,70},{82,90}},  rotation=0)));
  MixingVolumes.MixingVolume vol1(
    redeclare package Medium = Medium,
    V=100,
    m_flow_nominal=1,
    nPorts=2,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial) "Mixing volume"
                          annotation (Placement(transformation(extent={{100,80},
            {120,100}},rotation=0)));
  Sources.TraceSubstancesFlowSource sou1(
                                      redeclare package Medium = Medium,
      use_m_flow_in=true)
    annotation (Placement(transformation(extent={{-46,70},{-26,90}},  rotation=
            0)));
  Buildings.Utilities.Diagnostics.AssertEquality assEqu(threShold=1E-4)
    "Assert that both volumes have the same concentration"
    annotation (Placement(transformation(extent={{210,128},{230,148}},
                                                                     rotation=0)));
  MixingVolumes.MixingVolume vol2(
    redeclare package Medium = Medium,
    p_start=Medium.p_default,
    V=100,
    m_flow_nominal=1,
    nPorts=3,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial) "Mixing volume"
                          annotation (Placement(transformation(extent={{90,-20},
            {110,0}},  rotation=0)));
  MixingVolumes.MixingVolume vol3(
    redeclare package Medium = Medium,
    p_start=Medium.p_default,
    V=100,
    m_flow_nominal=1,
    nPorts=3,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial) "Mixing volume"
                          annotation (Placement(transformation(extent={{88,-60},
            {108,-40}},rotation=0)));
  Buildings.Utilities.Diagnostics.AssertEquality assEqu1(
                                                     threShold=1E-4)
    "Assert that both volumes have the same concentration"
    annotation (Placement(transformation(extent={{210,0},{230,20}},rotation=0)));
  MixingVolumes.MixingVolume vol4(
    redeclare package Medium = Medium,
    nPorts=4,
    p_start=Medium.p_default,
    V=100,
    m_flow_nominal=1,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial) "Mixing volume"
                          annotation (Placement(transformation(extent={{-16,-40},
            {4,-20}},  rotation=0)));
  Sources.TraceSubstancesFlowSource sou2(
                                      redeclare package Medium = Medium,
      use_m_flow_in=true)
    annotation (Placement(transformation(extent={{-48,-50},{-28,-30}}, rotation=
           0)));
  Buildings.Fluid.Sources.Boundary_pT bou(
    redeclare package Medium = Medium,
    p=101325,
    nPorts=1,
    T=293.15) annotation (Placement(transformation(extent={{-62,-80},{-42,-60}},
          rotation=0)));
  Buildings.Fluid.Sources.Boundary_pT sin(
    redeclare package Medium = Medium,
    nPorts=2,
    p=101320,
    T=293.15) "Sink boundary conditions"
              annotation (Placement(transformation(extent={{188,-50},{168,-30}},
          rotation=0)));
  FixedResistances.FixedResistanceDpM res1(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=1)
    "Resistance, used to check if species are transported between ports"
    annotation (Placement(transformation(extent={{126,-30},{148,-10}},rotation=
            0)));
  FixedResistances.FixedResistanceDpM res2(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=1)
    "Resistance, used to check if species are transported between ports"
    annotation (Placement(transformation(extent={{126,-70},{148,-50}},rotation=
            0)));
  FixedResistances.FixedResistanceDpM res3(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=1)
    "Resistance, used to check if species are transported between ports"
    annotation (Placement(transformation(extent={{-28,-80},{-6,-60}},  rotation=
           0)));
  inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial)
    annotation (Placement(transformation(extent={{-80,140},{-60,160}})));
  Sensors.TraceSubstances C(redeclare package Medium = Medium)
    "Trace substance sensor"
    annotation (Placement(transformation(extent={{120,134},{140,154}})));
  Sensors.TraceSubstances C1(redeclare package Medium = Medium)
    "Trace substance sensor"
    annotation (Placement(transformation(extent={{130,86},{150,106}})));
  Sensors.TraceSubstances C2(redeclare package Medium = Medium)
    "Trace substance sensor"
    annotation (Placement(transformation(extent={{168,6},{188,26}})));
  Sensors.TraceSubstances C3(redeclare package Medium = Medium)
    "Trace substance sensor"
    annotation (Placement(transformation(extent={{188,-50},{208,-30}})));
  FixedResistances.FixedResistanceDpM res4(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=1)
    "Resistance, used to check if species are transported between ports"
    annotation (Placement(transformation(extent={{58,-30},{80,-10}},  rotation=
            0)));
  FixedResistances.FixedResistanceDpM res6(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=1)
    "Resistance, used to check if species are transported between ports"
    annotation (Placement(transformation(extent={{58,-70},{80,-50}},  rotation=
            0)));
equation
  connect(res3.port_b, vol4.ports[2])
                                     annotation (Line(points={{-6,-70},{-6,-40},
          {-7,-40}},                                  color={0,127,255}));
  connect(res1.port_b, sin.ports[1])  annotation (Line(
      points={{148,-20},{158,-20},{158,-38},{168,-38}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(res2.port_b, sin.ports[2])  annotation (Line(
      points={{148,-60},{158,-60},{158,-42},{168,-42}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(bou.ports[1], res3.port_a) annotation (Line(
      points={{-42,-70},{-28,-70}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sou1.ports[1], res.port_a) annotation (Line(
      points={{-26,80},{60,80}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sou2.ports[1], vol4.ports[1]) annotation (Line(
      points={{-28,-40},{-9,-40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(step.y, sou.m_flow_in) annotation (Line(
      points={{-71,40},{-60,40},{-60,108},{-48.1,108}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, sou1.m_flow_in) annotation (Line(
      points={{-71,40},{-60,40},{-60,80},{-48.1,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(step.y, sou2.m_flow_in) annotation (Line(
      points={{-71,40},{-60,40},{-60,-40},{-50.1,-40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(assEqu.u1, C.C) annotation (Line(
      points={{208,144},{141,144}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(C1.C, assEqu.u2) annotation (Line(
      points={{151,96},{179.5,96},{179.5,132},{208,132}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(assEqu1.u1, C2.C) annotation (Line(
      points={{208,16},{189,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(C3.C, assEqu1.u2) annotation (Line(
      points={{209,-40},{220,-40},{220,-20},{200,-20},{200,4},{208,4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sou.ports[1], vol.ports[1]) annotation (Line(
      points={{-26,108},{108,108}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(vol.ports[2], C.port) annotation (Line(
      points={{112,108},{130,108},{130,134}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(res.port_b, vol1.ports[1]) annotation (Line(
      points={{82,80},{108,80}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(vol1.ports[2], C1.port) annotation (Line(
      points={{112,80},{140,80},{140,86}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(vol2.ports[1], res1.port_a) annotation (Line(
      points={{97.3333,-20},{126,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(vol3.ports[1], res2.port_a) annotation (Line(
      points={{95.3333,-60},{126,-60}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(C2.port, vol2.ports[2]) annotation (Line(
      points={{178,6},{178,0},{118,0},{118,-20},{100,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(C3.port, vol3.ports[2]) annotation (Line(
      points={{198,-50},{198,-78},{98,-78},{98,-60}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(vol4.ports[3], res4.port_a) annotation (Line(
      points={{-5,-40},{26,-40},{26,-20},{58,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(vol4.ports[4], res6.port_a) annotation (Line(
      points={{-3,-40},{26,-40},{26,-60},{58,-60}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(res6.port_b, vol3.ports[3]) annotation (Line(
      points={{80,-60},{100.667,-60}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(res4.port_b, vol2.ports[3]) annotation (Line(
      points={{80,-20},{92,-20},{92,-20},{102.667,-20}},
      color={0,127,255},
      smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,
            -100},{240,180}}), graphics),
            experiment(StopTime=600),
             __Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Fluid/Sources/Examples/TraceSubstancesFlowSource.mos"
        "Simulate and plot"),
    Documentation(info="<html>
This model demonstrates the use of trace substances that are added
to a volume of air.
The source is a step function of <i>2</i> kg/s CO<sub>2</sub> from <i>t=0</i> second
to <i>t=0.5</i> second.
The sensors <code>C</code> and <code>C1</code> measure the same concentration that initially increases
and then remains constant as there is no flow through the volumes <code>vol</code> and <code>vol1</code>.
The sensors 
<code>C2</code> and
<code>C3</code> first meaure an increase in concentration, which then decays to zero
as there is a mass flow rate with zero CO<sub>2</sub> from the source <code>bou</code> to the sink <code>sin</code>.
</html>", revisions="<html>
<ul>
<li>
September 19, 2013, by Michael Wetter:<br/>
Simplified example.
</li>
<li>
April 29, 2013, by Michael Wetter:<br/>
Changed the initialization of the medium volumes from free initial conditions
to <code>Modelica.Fluid.Types.Dynamics.FixedInitial</code>.
This was required for <code>vol</code> and <code>vol1</code> to have the same
computations for the initial states.
</li>
<li>
September 18, 2008 by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
end TraceSubstancesFlowSource;
