within Buildings.HeatTransfer.Interfaces;
package HAM
  partial package Interfaces

    partial model HAMDipole "partial model for heat and mass conductor"

    parameter Modelica.SIunits.Area A "Heat and mass transfer area";
      final parameter Modelica.SIunits.CoefficientOfHeatTransfer U = UA/A
        "U-value (without surface heat transfer coefficients)";
      final parameter Modelica.SIunits.ThermalConductance UA = 1/R
        "Thermal conductance of construction (without surface heat transfer coefficients)";
      parameter Modelica.SIunits.ThermalResistance R
        "Thermal resistance of construction";

        final parameter Real G = GA/A "G-Value";

        final parameter Real GA = 1/Rv "Vapour conductance of construction ";

        parameter Real Rv "Vapour propagation resistance of construction";

    Modelica.SIunits.TemperatureDifference dT;
    Modelica.SIunits.PressureDifference dp;

    public
      HAMPort_a port_a
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
      HAMPort_b port_b
        annotation (Placement(transformation(extent={{50,-10},{70,10}})));

    equation
      dT = port_a.T - port_b.T;
      dp = port_a.ham.p - port_b.ham.p;

      connect(port_a, port_b) annotation (Line(
          points={{-60,0},{60,0}},
          color={127,0,0},
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics));
    end HAMDipole;

    connector MassPort

        Modelica.SIunits.MassFraction Xi_outflow
        "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";

      flow Modelica.SIunits.MassFlowRate m_flow
        "Mass flow rate from the connection point into the component";

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Rectangle(
              extent={{-100,44},{100,-80}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid)}), uses(Modelica(version="3.2.1")),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={Rectangle(extent={{-100,100},{100,-100}},
                lineColor={0,127,0},
              fillPattern=FillPattern.Solid,
              fillColor={0,0,255})}));

    end MassPort;

    connector HeatMassPort_a

      parameter Boolean use_massPort = true;

      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort;

      MassPort massPort if use_massPort;

      annotation (uses(Modelica(version="3.2.1")), Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Polygon(
              points={{0,100},{-72,40},{0,100}},
              pattern=LinePattern.None,
              smooth=Smooth.None,
              fillColor={98,0,98},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{0,100},{0,100}},
              pattern=LinePattern.None,
              smooth=Smooth.None,
              fillColor={98,0,98},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{0,100},{-98,-2},{0,100}},
              pattern=LinePattern.None,
              smooth=Smooth.None,
              fillColor={98,0,98},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{0,100},{0,100}},
              pattern=LinePattern.None,
              smooth=Smooth.None,
              fillColor={98,0,98},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{0,100},{-98,0},{2,-100},{100,0},{0,100}},
              pattern=LinePattern.None,
              smooth=Smooth.None,
              fillColor={98,0,98},
              fillPattern=FillPattern.Solid)}));
    end HeatMassPort_a;

    connector HeatMassPort_b

      parameter Boolean use_massPort = true;

      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort;

      MassPort massPort if use_massPort;

      annotation (uses(Modelica(version="3.2.1")), Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{-2,100},{-100,0},{0,-100},{98,0},{-2,100}},
              pattern=LinePattern.None,
              smooth=Smooth.None,
              fillColor={98,0,98},
              fillPattern=FillPattern.Solid)}));
    end HeatMassPort_b;
  end Interfaces;

  package Composants

    model HAMConductor
      extends Interfaces.HAMDipole( final R = material.x/material.lamda/A, Rv = material.x/material.delta/A);

    Modelica.SIunits.Pressure p[nSta];
    Modelica.SIunits.MassFlowRate m_flow[nSta+1];
    Modelica.SIunits.HeatFlowRate Q_flow[nSta+1];
    Modelica.SIunits.Temperature T[nSta];

      replaceable parameter
        Buildings.HeatTransfer.Interfaces.HAM.Media.GenericHAM                     material
        annotation (Evaluate=true, choicesAllMatching=true, Placement(transformation(extent={{60,60},
                {80,80}})));

    protected
      final parameter Integer nSta(min=1) = material.nSta
        "Number of state variables";
      final parameter Modelica.SIunits.ThermalConductance UAnSta = UA*nSta
        "Thermal conductance between nodes";
      final parameter Real GAnSta = GA*nSta;

      final parameter Modelica.SIunits.ThermalConductance UAnSta2 = 2*UAnSta
        "Thermal conductance between nodes and surface boundary";
      final parameter Real GAnSta2 = 2*GAnSta;

      parameter Modelica.SIunits.Mass m = A*material.x*material.rho/material.nSta
        "Mass associated with the temperature state";
      parameter Modelica.SIunits.HeatCapacity C = m*material.cv
        "Heat capacity associated with the temperature state";

    equation
      port_a.Q_flow = +Q_flow[1];
      port_b.Q_flow = -Q_flow[nSta+1];

      port_a.ham.m_flow = +m_flow[1];
      port_b.ham.m_flow = -m_flow[nSta+1];

        port_a.T-T[1] = Q_flow[1]/UAnSta2;
        T[nSta] -port_b.T = Q_flow[nSta+1]/UAnSta2;

        port_a.ham.p-p[1] = m_flow[1]/GAnSta2;
        p[nSta] -port_b.ham.p = m_flow[nSta+1]/GAnSta2;

        for i in 2:nSta loop
           // Q_flow[i] is heat flowing from (i-1) to (i)
           T[i-1]-T[i] = Q_flow[i]/UAnSta;
           p[i-1]-p[i] = m_flow[i]/GAnSta;

        end for;

        for i in 1:nSta loop
          der(T[i]) = (Q_flow[i]-Q_flow[i+1])/C;
          der(p[i]) = (m_flow[i]-m_flow[i+1]);

        end for;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}}), graphics={
            Rectangle(
              extent={{-40,84},{-28,-86}},
              lineColor={0,0,0},
              fillColor={175,175,175},
              fillPattern=FillPattern.CrossDiag),
            Rectangle(
              extent={{28,84},{40,-86}},
              lineColor={0,0,0},
              fillColor={175,175,175},
              fillPattern=FillPattern.CrossDiag),
            Polygon(
              points={{12,10},{14,10},{16,6},{18,0},{18,-4},{16,-10},{10,-14},{
                  6,-18},{-2,-20},{-6,-16},{-12,-10},{-14,0},{-12,6},{-10,10},{
                  -8,12},{-6,14},{-2,16},{2,16},{8,14},{12,10}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-52,4},{52,-4}},
              lineColor={0,0,127},
              fillColor={175,8,19},
              fillPattern=FillPattern.HorizontalCylinder),
            Line(
              points={{-28,60},{30,60}},
              color={0,0,255},
              pattern=LinePattern.Dot,
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.Filled}),
            Line(
              points={{-28,-60},{28,-60},{32,-60}},
              color={0,0,255},
              pattern=LinePattern.Dot,
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.Filled})}));
    end HAMConductor;

    model HAMConductor3

    Modelica.SIunits.MassFlowRate m_flow[nSta+1];
    Modelica.SIunits.HeatFlowRate Q_flow[nSta+1];
    Modelica.SIunits.Temperature T[nSta](start=
         {T_a_start+(T_b_start-T_a_start) * UA *
            sum(1/(if (k==1 or k==nSta+1) then UAnSta2 else UAnSta) for k in 1:i) for i in 1:nSta},
          each nominal = 300) "Temperature at the states";

    Real D[nSta]
        "water vapour diffusion coefficient in air (Unit = Kg.m^-2.s^-1)";
    Modelica.SIunits.Pressure pw[nSta] "Vapour Pressure at the state";
    Real phi[nSta](each start=0.5) "relative humidity";
    Real dw_dphi[nSta];

    Modelica.SIunits.Pressure psat[nSta]
        "equilibrum vapor pressure at the states";
    Modelica.SIunits.MassFraction Xi_outflow[nSta] "mass fraction at the state";

      replaceable parameter Data.BaseClasses.HygroThermalMaterial material
        "Material from Data.Solids, Data.SolidsPCM or Data.Resistances"
        annotation (Evaluate=true, choicesAllMatching=true, Placement(transformation(extent={{60,60},
                {80,80}})));

    protected
      parameter Modelica.SIunits.Area A= 1 "Heat transfer area";
      final parameter Modelica.SIunits.CoefficientOfHeatTransfer U = UA/A
        "U-value (without surface heat transfer coefficients)";
      parameter Modelica.SIunits.ThermalResistance R = material.x/(material.k*A)
        "Thermal resistance of construction";
      final parameter Modelica.SIunits.ThermalConductance UA = 1/R
        "Thermal conductance of construction (without surface heat transfer coefficients)";

      Modelica.SIunits.TemperatureDifference dT "port_a.T - port_b.T";

      final parameter Integer nSta(min=1) = material.nSta
        "Number of state variables";
      final parameter Modelica.SIunits.ThermalConductance UAnSta = UA*nSta
        "Thermal conductance between nodes";
      final parameter Modelica.SIunits.ThermalConductance UAnSta2 = 2*UAnSta
        "Thermal conductance between nodes and surface boundary";
      parameter Modelica.SIunits.Mass m = A*material.x*material.d/material.nSta
        "Mass associated with the temperature state";
      parameter Modelica.SIunits.HeatCapacity C = m*material.c
        "Heat capacity associated with the temperature state";
      parameter Real  b = material.b;
      parameter Modelica.SIunits.Pressure patm = 101325 "atmospheric pressure";

      parameter Modelica.SIunits.Temperature T_a_start=293.15
        annotation (Dialog(group="Initialization", enable=not steadyStateInitial));

      parameter Modelica.SIunits.Temperature T_b_start=293.15
        annotation (Dialog(group="Initialization", enable=not steadyStateInitial));

    public
      Interfaces.HeatMassPort_b heatMassPort_b
        annotation (Placement(transformation(extent={{44,-10},{64,10}}),
            iconTransformation(extent={{44,-10},{64,10}})));
      Interfaces.HeatMassPort_a heatMassPort_a annotation (Placement(transformation(
              extent={{-68,-10},{-48,10}}), iconTransformation(extent={{-68,-10},{-48,
                10}})));
    equation

      dT = heatMassPort_a.heatPort.T - heatMassPort_b.heatPort.T;

      psat[1]= Modelica.Math.exp(-5800/heatMassPort_a.heatPort.T+1.391-0.04864*heatMassPort_a.heatPort.T+4.176*10^(-5)*heatMassPort_a.heatPort.T^2-1.445*10^(-8)*heatMassPort_a.heatPort.T^3+6.545*Modelica.Math.log(heatMassPort_a.heatPort.T));
      //phi[1]= (heatMassPort_a.massPort.Xi_outflow*patm)/(heatMassPort_a.massPort.Xi_outflow*psat[1]+0.62198*psat[1]);

      //Xi_outflow[1]= heatMassPort_a.massPort.Xi_outflow;

      heatMassPort_a.massPort.m_flow = +m_flow[1];
      Buildings.Utilities.Psychrometrics.Functions.pW_X(heatMassPort_a.massPort.Xi_outflow) - pw[1]= m_flow[1]/(nSta*D[1]*A/(material.x*material.Mu));

      heatMassPort_a.heatPort.Q_flow =+Q_flow[1];
      heatMassPort_a.heatPort.T-T[1] = Q_flow[1]/UAnSta2;

      heatMassPort_b.massPort.m_flow= -m_flow[nSta+1];
      pw[nSta] - Buildings.Utilities.Psychrometrics.Functions.pW_X(heatMassPort_b.massPort.Xi_outflow) = m_flow[nSta+1]/(nSta*D[nSta]*A/(material.x*material.Mu));

      heatMassPort_b.heatPort.Q_flow = -Q_flow[nSta+1];
      T[nSta] - heatMassPort_b.heatPort.T = Q_flow[nSta+1]/UAnSta2;

      for i in 1:nSta loop
      pw[i]=psat[i]*phi[i];
        D[i]= (2*10^(-7)*T[i]^0.81)/patm;
        dw_dphi[i]= material.xw_f*(b-1)*(b+1-2*phi[i])/((b-phi[i])*(b-phi[i]));

     // phi[i]= (Xi_outflow[i]*patm)/(0.62198*psat[i]+Xi_outflow[i]*psat[i]);
      Xi_outflow[i] =Buildings.Utilities.Psychrometrics.Functions.X_pW(pw[i]);
    end for;

      for i in 2:nSta loop

        psat[i]= Modelica.Math.exp(-5800/T[i]+1.391-0.04864*T[i]+4.176*10^(-5)*T[i]^2-1.445*10^(-8)*T[i]^3+6.545*Modelica.Math.log(T[i]));
        T[i-1]-T[i] = Q_flow[i]/UAnSta;
        pw[i-1]-pw[i] = m_flow[i]/((nSta*D[i]*A)/(material.x*material.Mu));

      end for;

      for i in 1:nSta loop
          der(T[i]) = (Q_flow[i]-Q_flow[i+1])/C;
          der(phi[i]) = ((m_flow[i]-m_flow[i+1])/dw_dphi[i])/(A/(0.4/nSta));
      end for;

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics), Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              extent={{-40,84},{-28,-86}},
              lineColor={0,0,0},
              fillColor={175,175,175},
              fillPattern=FillPattern.CrossDiag),
            Rectangle(
              extent={{28,82},{40,-88}},
              lineColor={0,0,0},
              fillColor={175,175,175},
              fillPattern=FillPattern.CrossDiag),
            Rectangle(
              extent={{-54,4},{50,-4}},
              lineColor={0,0,127},
              fillColor={175,8,19},
              fillPattern=FillPattern.HorizontalCylinder),
            Polygon(
              points={{10,12},{12,12},{14,8},{16,2},{16,-2},{14,-8},{8,-12},{4,-16},
                  {-4,-18},{-8,-14},{-14,-8},{-16,2},{-14,8},{-12,12},{-10,14},{-8,16},
                  {-4,18},{0,18},{6,16},{10,12}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-30,60},{28,60}},
              color={0,0,255},
              pattern=LinePattern.Dot,
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.Filled}),
            Line(
              points={{-30,-60},{26,-60},{30,-60}},
              color={0,0,255},
              pattern=LinePattern.Dot,
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.Filled})}));
    end HAMConductor3;

    model PrescribedHumidity

      Modelica.Blocks.Interfaces.RealInput X annotation (Placement(transformation(
              extent={{-140,-20},{-100,20}}, rotation=0)));
      Interfaces.MassPort massPort annotation (Placement(transformation(extent={{96,
                -10},{116,10}}), iconTransformation(extent={{96,-10},{116,10}})));
    equation
      massPort.Xi_outflow = X;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
                             graphics={Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,127},
              fillColor={170,170,255},
              fillPattern=FillPattern.Backward), Line(
              points={{-100,0},{100,0},{102,0}},
              color={0,0,127},
              pattern=LinePattern.Dash,
              thickness=1,
              smooth=Smooth.None,
              arrow={Arrow.None,Arrow.Filled})}));
    end PrescribedHumidity;

    model FixedHumidity
     parameter Modelica.SIunits.MassFraction X "Fixed humidity at port";

      Interfaces.MassPort massPort
        annotation (Placement(transformation(extent={{80,-8},{100,12}})));

    equation
      massPort.Xi_outflow = X;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}}), graphics={Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,127},
              fillColor={170,170,255},
              fillPattern=FillPattern.Backward), Line(
              points={{-60,4},{90,2}},
              color={0,0,127},
              pattern=LinePattern.Dash,
              smooth=Smooth.None,
              arrow={Arrow.None,Arrow.Filled},
              thickness=1)}));
    end FixedHumidity;

    model HAMConductor4

    Modelica.SIunits.MassFlowRate m_flow[nSta+1];
    Modelica.SIunits.HeatFlowRate Q_flow[nSta+1];
    Modelica.SIunits.Temperature T[nSta](start=
         {T_a_start+(T_b_start-T_a_start) * UA *
            sum(1/(if (k==1 or k==nSta+1) then UAnSta2 else UAnSta) for k in 1:i) for i in 1:nSta},
          each nominal = 300) "Temperature at the states";

    Real D[nSta]
        "water vapour diffusion coefficient in air (Unit = Kg.m^-2.s^-1)";
    Modelica.SIunits.Pressure pw[nSta] "water vapour pressure";
    Real phi[nSta](each start=0.5) "relative humidity";
    Real dw_dphi[nSta];
    Modelica.SIunits.Pressure psat[nSta]
        "equilibrum vapor pressure at the states";
    Modelica.SIunits.MassFraction Xi_outflow[nSta] "mass fraction at the state";
    Real w[nSta] "water content at the states";

      replaceable parameter Data.BaseClasses.HygroThermalMaterial material
        "Material from Data.Solids, Data.SolidsPCM or Data.Resistances"
        annotation (Evaluate=true, choicesAllMatching=true, Placement(transformation(extent={{60,60},
                {80,80}})));

    protected
      parameter Modelica.SIunits.Area A= 1 "Heat transfer area";
      final parameter Modelica.SIunits.CoefficientOfHeatTransfer U = UA/A
        "U-value (without surface heat transfer coefficients)";
      parameter Modelica.SIunits.ThermalResistance R = material.x/(material.k*A)
        "Thermal resistance of construction";
      final parameter Modelica.SIunits.ThermalConductance UA = 1/R
        "Thermal conductance of construction (without surface heat transfer coefficients)";

      Modelica.SIunits.TemperatureDifference dT "port_a.T - port_b.T";

      final parameter Integer nSta(min=1) = material.nSta
        "Number of state variables";
      final parameter Modelica.SIunits.ThermalConductance UAnSta = UA*nSta
        "Thermal conductance between nodes";
      final parameter Modelica.SIunits.ThermalConductance UAnSta2 = 2*UAnSta
        "Thermal conductance between nodes and surface boundary";
      parameter Modelica.SIunits.Mass m = A*material.x*material.d/material.nSta
        "Mass associated with the temperature state";
      parameter Modelica.SIunits.HeatCapacity C = m*material.c
        "Heat capacity associated with the temperature state";
      parameter Real  b = material.b;
      parameter Modelica.SIunits.Pressure patm = 101325 "atmospheric pressure";

      parameter Modelica.SIunits.Temperature T_a_start=293.15
        "Initial temperature at port_a, used if steadyStateInitial = false"
        annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
      parameter Modelica.SIunits.Temperature T_b_start=293.15
        "Initial temperature at port_b, used if steadyStateInitial = false"
        annotation (Dialog(group="Initialization", enable=not steadyStateInitial));

    public
      Interfaces.HeatMassPort_b heatMassPort_b
        annotation (Placement(transformation(extent={{44,-10},{64,10}}),
            iconTransformation(extent={{44,-10},{64,10}})));
      Interfaces.HeatMassPort_a heatMassPort_a annotation (Placement(transformation(
              extent={{-68,-10},{-48,10}}), iconTransformation(extent={{-68,-10},{-48,
                10}})));
    equation

      dT = heatMassPort_a.heatPort.T - heatMassPort_b.heatPort.T;

      psat[1]= Modelica.Math.exp(-5800/heatMassPort_a.heatPort.T+1.391-0.04864*heatMassPort_a.heatPort.T+4.176*10^(-5)*heatMassPort_a.heatPort.T^2-1.445*10^(-8)*heatMassPort_a.heatPort.T^3+6.545*Modelica.Math.log(heatMassPort_a.heatPort.T));
      phi[1]= (heatMassPort_a.massPort.Xi_outflow*patm)/(heatMassPort_a.massPort.Xi_outflow*psat[1]+0.62198*psat[1]);

      Xi_outflow[1]= heatMassPort_a.massPort.Xi_outflow;
      Xi_outflow[nSta] = heatMassPort_a.massPort.Xi_outflow;

      heatMassPort_a.heatPort.Q_flow =+Q_flow[1];
      heatMassPort_a.heatPort.T-T[1] = Q_flow[1]/UAnSta2;

      heatMassPort_a.massPort.m_flow = +m_flow[1];

      heatMassPort_b.heatPort.Q_flow = -Q_flow[nSta+1];
      T[nSta] - heatMassPort_b.heatPort.T = Q_flow[nSta+1]/UAnSta2;

      heatMassPort_b.massPort.m_flow = -m_flow[nSta+1];

      for i in 1:nSta loop
        pw[i]=phi[i]*psat[i];
        D[i]= (2*10^(-7)*T[i]^0.81)/patm;
        dw_dphi[i]= material.xw_f*(b-1)*(b+1-2*phi[i])/((b-phi[i])*(b-phi[i]));
        w[i]=material.xw_f*((b-1)*phi[i])/(b-phi[i]);
      end for;

      for i in 2:nSta loop

        phi[i]= (Xi_outflow[i]*patm)/(0.62198*psat[i]+Xi_outflow[i]*psat[i]);
        psat[i]= Modelica.Math.exp(-5800/T[i]+1.391-0.04864*T[i]+4.176*10^(-5)*T[i]^2-1.445*10^(-8)*T[i]^3+6.545*Modelica.Math.log(T[i]));
        T[i-1]-T[i] = Q_flow[i]/UAnSta;
        phi[i-1]*psat[i-1]-phi[i]*psat[i] = m_flow[i]/(nSta*D[i]*A/(material.x*material.Mu));

      end for;

      for i in 1:nSta loop
          der(T[i]) = (Q_flow[i]-Q_flow[i+1])/C;
          der(phi[i]) = (m_flow[i]-m_flow[i+1])/dw_dphi[i];
      end for;

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics), Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              extent={{-40,84},{-28,-86}},
              lineColor={0,0,0},
              fillColor={175,175,175},
              fillPattern=FillPattern.CrossDiag),
            Rectangle(
              extent={{28,82},{40,-88}},
              lineColor={0,0,0},
              fillColor={175,175,175},
              fillPattern=FillPattern.CrossDiag),
            Rectangle(
              extent={{-54,4},{50,-4}},
              lineColor={0,0,127},
              fillColor={175,8,19},
              fillPattern=FillPattern.HorizontalCylinder),
            Polygon(
              points={{10,12},{12,12},{14,8},{16,2},{16,-2},{14,-8},{8,-12},{4,-16},
                  {-4,-18},{-8,-14},{-14,-8},{-16,2},{-14,8},{-12,12},{-10,14},{-8,16},
                  {-4,18},{0,18},{6,16},{10,12}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-30,60},{28,60}},
              color={0,0,255},
              pattern=LinePattern.Dot,
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.Filled}),
            Line(
              points={{-30,-60},{26,-60},{30,-60}},
              color={0,0,255},
              pattern=LinePattern.Dot,
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.Filled})}));
    end HAMConductor4;
  end Composants;

  package Media "medium properties"
    extends Modelica.Icons.MaterialPropertiesPackage;

    record GenericHAM "Record containing HAM material properties"
      extends Modelica.Icons.Record;
      parameter Modelica.SIunits.Length x = 0.2 "Material thickness";
      parameter Modelica.SIunits.Density rho = 1 "Density";
      parameter Modelica.SIunits.SpecificHeatCapacity cp = 1
        "Specific heat capacity at constant pressure";
      parameter Modelica.SIunits.SpecificHeatCapacity cv = 1
        "Specific heat capacity at constant volume";
      parameter Modelica.SIunits.ThermalConductivity lamda = 1
        "Thermal conductivity";
      parameter Real R=1 "Thermal resistance of a unit area of material";
      parameter Integer nSta = 3 "Actual number of state variables in material";
      parameter Real delta = 1 "Vapour permeability";
      annotation (Documentation(info="<html>
Record containing material properties.
</html>"));
    end GenericHAM;

    record ConcreteHygroThermal =
        Buildings.HeatTransfer.Interfaces.HAM.Media.GenericHAM (
        rho=2230,
        cp=2240,
        cv=840,
        lamda = 1.6,
        R = 0,
        nSta= 5,
        delta = 0.00005) "Concrete ";
  end Media;

  model test
  //package Medium = Buildings.Media.GasesConstantDensity.MoistAirUnsaturated
    //  "Air model used in the example model";
    Composants.HAMConductor3 hAMConductor3_1(redeclare
        Buildings.HeatTransfer.Data.Solids.ConcreteHAM material(x=0.4))
      annotation (Placement(transformation(extent={{0,-30},{50,22}})));
    Interfaces.HeatMassPort_a heatMassPort_a1
      annotation (Placement(transformation(extent={{-26,-6},{-14,6}})));
    Interfaces.HeatMassPort_b heatMassPort_b1
      annotation (Placement(transformation(extent={{54,-4},{62,4}})));
    Sources.FixedTemperature fixedTemperature(T=283.15)
      annotation (Placement(transformation(extent={{-86,-12},{-66,8}})));
    Sources.FixedTemperature fixedTemperature1(T=293.15) annotation (Placement(
          transformation(
          extent={{-10,10},{10,-10}},
          rotation=180,
          origin={84,8})));
    Composants.FixedHumidity fixedHumidity1(X=0.01)
      annotation (Placement(transformation(extent={{30,50},{50,70}})));
    Composants.FixedHumidity fixedHumidity(X=0.004)
      annotation (Placement(transformation(extent={{-66,26},{-46,46}})));
  equation
    //hAMConductor3_1.heatMassPort_a.massPort.m_flow = prescribedHumidity.massPort.m_flow;
      //  connect( prescribedTemperature.port, hAMConductor3_1.heatMassPort_a.heatPort);
      //  connect( prescribedHumidity.massPort, hAMConductor3_1.heatMassPort_a.massPort);

    connect(hAMConductor3_1.heatMassPort_a, heatMassPort_a1) annotation (Line(
        points={{10.5,-4},{-6,-4},{-6,0},{-20,0}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(hAMConductor3_1.heatMassPort_b, heatMassPort_b1) annotation (Line(
        points={{38.5,-4},{48,-4},{48,0},{58,0}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(heatMassPort_a1, heatMassPort_a1) annotation (Line(
        points={{-20,0},{-20,0}},
        color={0,0,0},
        pattern=LinePattern.None,
        smooth=Smooth.None));
    connect(fixedTemperature.port, heatMassPort_a1.heatPort) annotation (Line(
        points={{-66,-2},{-44,-2},{-44,0},{-20,0}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedTemperature1.port, heatMassPort_b1.heatPort) annotation (Line(
        points={{74,8},{66,8},{66,0},{58,0}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedHumidity1.massPort, heatMassPort_b1.massPort) annotation (Line(
        points={{49,60.2},{49,31.1},{58,31.1},{58,0}},
        color={0,127,0},
        smooth=Smooth.None));
    connect(fixedHumidity.massPort, heatMassPort_a1.massPort) annotation (Line(
        points={{-47,36.2},{-47,16.1},{-20,16.1},{-20,0}},
        color={0,127,0},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),      graphics));
  end test;

  model test2
    Composants.HAMConductor4 hAMConductor4_1(redeclare
        Buildings.HeatTransfer.Data.Solids.ConcreteHAM material(x=0.4), phi(
          fixed=true))
      annotation (Placement(transformation(extent={{-22,-24},{22,24}})));
    Interfaces.HeatMassPort_a heatMassPort_a1
      annotation (Placement(transformation(extent={{-46,-6},{-34,6}})));
    Interfaces.HeatMassPort_b heatMassPort_b1
      annotation (Placement(transformation(extent={{34,-6},{46,6}})));
    Composants.FixedHumidity fixedHumidity(X=0.002)
      annotation (Placement(transformation(extent={{-88,-34},{-68,-14}})));
    Sources.FixedTemperature fixedTemperature(T=283.15)
      annotation (Placement(transformation(extent={{-90,22},{-70,42}})));
    Sources.FixedTemperature fixedTemperature1(T=293.15)
      annotation (Placement(transformation(extent={{74,22},{54,42}})));
  equation
    connect(hAMConductor4_1.heatMassPort_a, heatMassPort_a1) annotation (Line(
        points={{-12.76,0},{-40,0}},
        color={0,0,0},
        pattern=LinePattern.None,
        smooth=Smooth.None));
    connect(hAMConductor4_1.heatMassPort_b, heatMassPort_b1) annotation (Line(
        points={{11.88,0},{40,0}},
        color={0,0,0},
        pattern=LinePattern.None,
        smooth=Smooth.None));
    connect(fixedHumidity.massPort, heatMassPort_a1.massPort) annotation (Line(
        points={{-69,-23.8},{-54.5,-23.8},{-54.5,0},{-40,0}},
        color={0,127,0},
        smooth=Smooth.None));
    connect(fixedTemperature.port, heatMassPort_a1.heatPort) annotation (Line(
        points={{-70,32},{-54,32},{-54,0},{-40,0}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(fixedTemperature1.port, heatMassPort_b1.heatPort) annotation (Line(
        points={{54,32},{48,32},{48,0},{40,0}},
        color={191,0,0},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),      graphics));
  end test2;
end HAM;
