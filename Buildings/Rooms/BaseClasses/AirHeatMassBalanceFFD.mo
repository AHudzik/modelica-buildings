within Buildings.Rooms.BaseClasses;
model AirHeatMassBalanceFFD
  "Heat and mass balance of the air based on fast fluid flow dynamics"
  extends Buildings.Rooms.BaseClasses.PartialAirHeatMassBalance(
   energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial,
   massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial);

  parameter Boolean useFFD = true
    "Set to false to deactivate the FFD interface and use instead yFixed as output"
    annotation(Evaluate = true);

  parameter Modelica.SIunits.Time samplePeriod(min=100*Modelica.Constants.eps)
    "Sample period of component"
    annotation(Dialog(group = "Sampling"));
  parameter Modelica.SIunits.Time startTime
    "First sample time instant. fixme: this should be at first step."
    annotation(Dialog(group = "Sampling"));

  parameter Boolean haveSensor
    "Flag, true if the model has at least one sensor";
  parameter Integer nSen(min=0)
    "Number of sensors that are connected to CFD output";
  parameter String sensorName[nSen]
    "Names of sensors as declared in the CFD input file";

  // fixme: for the ffd instance, need to correctly assign flaWri
  FFDExchange ffd(
    final startTime=startTime,
    final activateInterface=useFFD,
    final samplePeriod = if useFFD then samplePeriod else Modelica.Constants.inf,
    final uStart=uStart,
    final nWri=kFluIntC_inflow+Medium.nC*nPorts,
    final nRea=kSen+nSen,
    final nSur=nSur,
    final surIde = surIde,
    final haveShade = haveShade,
    final haveSensor=haveSensor,
    final nSen=nSen,
    final sensorName=sensorName,
    final yFixed=yFixed) "Block that exchanges data with the FFD simulation"
    annotation (Placement(transformation(extent={{-40,180},{-20,200}})));

  Modelica.Blocks.Interfaces.RealOutput yCFD[nSen] if
      haveSensor "Sensor for output from CFD"
    annotation (Placement(transformation(
     extent={{-10,-10},{10,10}},
        rotation=270,
        origin={180,-250})));

  // Values that are used for uStart
protected
  parameter Real uStart[kFluIntC_inflow+Medium.nC*nPorts](fixed=false)
    "Values used for uStart in FFDExchange";

  // Values that are used for yFixed
  parameter Real yFixed[kSen+nSen](fixed=false)
    "Values used for yFixed in FFDExchange";

  parameter Modelica.SIunits.HeatFlowRate Q_flow_fixed[kSurBou+nSurBou]=fill(0, kSurBou+nSurBou)
    "Surface heat flow rate used for yFixed"
    annotation (Dialog(group="Outputs if activateInterface=false"));
  parameter Modelica.SIunits.Temperature TRooAve_fixed = Medium.T_default
    "Average room air temperature used for yFixed"
    annotation (Dialog(group="Outputs if activateInterface=false"));
  parameter Modelica.SIunits.Temperature TSha_fixed[NConExtWin] = fill(T_start, NConExtWin)
    "Shade temperature used for yFixed"
    annotation (Dialog(group="Outputs if activateInterface=false"));
  parameter Modelica.SIunits.Temperature T_outflow_fixed[nPorts] = fill(T_start, nPorts)
    "Temperature of the fluid that flows into the HVAC system used for yFixed"
    annotation (Dialog(group="Outputs if activateInterface=false"));
  parameter Real Xi_outflow_fixed[nPorts*Medium.nXi](fixed=false)
    "Species concentration of the fluid that flows into the HVAC system used for yFixed"
    annotation (Dialog(group="Outputs if activateInterface=false"));
  parameter Real C_outflow_fixed[nPorts*max(1, Medium.nC)]=
    if Medium.nC == 0 then fill(0, nPorts) else cat(1, C_start[i] for i in 1:Medium.nC)
    "Trace substances of the fluid that flows into the HVAC system used for yFixed"
    annotation (Dialog(group="Outputs if activateInterface=false"));

  final parameter FFDSurfaceIdentifier surIde[kSurBou+nSurBou]=
   assignSurfaceIdentifier(
    nConExt=  nConExt,
    nConExtWin=  nConExtWin,
    nConPar=  nConPar,
    nConBou=  nConBou,
    nSurBou=  nSurBou,
    nSur=  nSur,
    haveShade=  haveShade,
    nameConExt=    datConExt.name,
    AConExt=       datConExt.A,
    tilConExt=     datConExt.til,
    bouConConExt=  datConExt.boundaryCondition,
    nameConExtWin=    datConExtWin.name,
    AConExtWin=       datConExtWin.AOpa,
    tilConExtWin=     datConExtWin.til,
    bouConConExtWin=  datConExtWin.boundaryCondition,
    AGla=       datConExtWin.AGla,
    AFra=       datConExtWin.AFra,
    nameConPar=    datConPar.name,
    AConPar=       datConPar.A,
    tilConPar=     datConPar.til,
    bouConConPar=  datConPar.boundaryCondition,
    nameConBou=    datConBou.name,
    AConBou=       datConBou.A,
    tilConBou=     datConBou.til,
    bouConConBou=  datConBou.boundaryCondition,
    nameSurBou=    surBou.name,
    ASurBou=       surBou.A,
    tilSurBou=     surBou.til,
    bouConSurBou=  surBou.boundaryCondition)
    "Names of all surfaces in the order in which their properties are sent to FFD"
     annotation (Evaluate=true);

  // Interfaces between the FFD block and the heat ports of this model
  // Here, we directly access datConExt instead of surIde. The reason is
  // the Dymola thinks that surIde.bouCon is not fixed at translation time
  // and then refuses to use this parameter to conditionally remove connectors
  // in FFDSurfaceInterface.
  FFDSurfaceInterface ffdConExt[NConExt](
    final bouCon={datConExt[i].boundaryCondition for i in 1:NConExt}) if
       haveConExt "Interface to heat port of exterior constructions"
    annotation (Placement(transformation(extent={{180,210},{200,230}})));

  FFDSurfaceInterface ffdConExtWin[NConExtWin](
    final bouCon={datConExtWin[i].boundaryCondition for i in 1:NConExtWin}) if
        haveConExtWin
    "Interface to heat port of opaque part of exterior constructions with window"
    annotation (Placement(transformation(extent={{180,170},{200,190}})));

  FFDSurfaceInterface ffdGlaUns[NConExtWin](
    final bouCon={datConExtWin[i].boundaryCondition for i in 1:NConExtWin}) if
       haveConExtWin "Interface to heat port of unshaded part of glass"
    annotation (Placement(transformation(extent={{180,110},{200,130}})));

  FFDSurfaceInterface ffdGlaSha[NConExtWin](
    final bouCon={datConExtWin[i].boundaryCondition for i in 1:NConExtWin}) if
       haveShade "Interface to heat port of shaded part of glass"
    annotation (Placement(transformation(extent={{180,70},{200,90}})));

  FFDSurfaceInterface ffdConExtWinFra[NConExtWin](
    final bouCon={datConExtWin[i].boundaryCondition for i in 1:NConExtWin}) if
       haveConExtWin "Interface to heat port of window frame"
    annotation (Placement(transformation(extent={{180,-10},{200,10}})));

  FFDSurfaceInterface ffdConPar_a[NConPar](
    final bouCon={datConPar[i].boundaryCondition for i in 1:NConPar}) if
       haveConPar
    "Interface to heat port of surface a of partition constructions"
    annotation (Placement(transformation(extent={{180,-70},{200,-50}})));

  FFDSurfaceInterface ffdConPar_b[NConPar](
    final bouCon={datConPar[i].boundaryCondition for i in 1:NConPar}) if
       haveConPar
    "Interface to heat port of surface b of partition constructions"
    annotation (Placement(transformation(extent={{180,-110},{200,-90}})));

  FFDSurfaceInterface ffdConBou[NConBou](
    final bouCon={datConBou[i].boundaryCondition for i in 1:NConBou}) if
       haveConBou
    "Interface to heat port that connects to room-side surface of constructions that expose their other surface to the outside"
    annotation (Placement(transformation(extent={{180,-170},{200,-150}})));

  FFDSurfaceInterface ffdSurBou[NSurBou](
    final bouCon={surBou[i].boundaryCondition for i in 1:NSurBou}) if
       haveSurBou
    "Interface to heat port of surfaces of models that compute the heat conduction outside of this room"
    annotation (Placement(transformation(extent={{180,-230},{200,-210}})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature ffdHeaPorAir
    "Interface to heat port of air node"
    annotation (Placement(transformation(extent={{-140,-10},{-160,10}})));

  FFDFluidInterface fluInt(
    redeclare final package Medium = Medium,
    final nPorts=nPorts,
    final V=V,
    final p_start=p_start,
    final T_start=T_start,
    final X_start=X_start,
    final C_start=C_start,
    final C_nominal=C_nominal) "Fluid interface"
    annotation (Placement(transformation(extent={{10,-198},{-10,-178}})));

  // The following list declares the first index minus 1
  // of the input and output signals to the FFD block.
  // These parameters are then used to loop over the connectors, such
  // as
  //    for i in kConExt+1:kConExt+nConExt loop
  //      ...
  //    end for;
  final parameter Integer kConExt = 0
    "Offset used to connect FFD signals to conExt";
  final parameter Integer kConExtWin = kConExt + nConExt
    "Offset used to connect FFD signals to conExtWin";
  final parameter Integer kGlaUns =    kConExtWin + nConExtWin
    "Offset used to connect FFD signals to glaUns";
  final parameter Integer kGlaSha = kGlaUns + nConExtWin
    "Offset used to connect FFD signals to glaSha";
  final parameter Integer kConExtWinFra = if haveShade then kGlaSha + nConExtWin else kGlaSha
    "Offset used to connect FFD signals to glaSha";
  final parameter Integer kConPar_a = kConExtWinFra + nConExtWin
    "Offset used to connect FFD signals to conPar_a";
  final parameter Integer kConPar_b = kConPar_a + nConPar
    "Offset used to connect FFD signals to conPar_b";
  final parameter Integer kConBou = kConPar_b + nConPar
    "Offset used to connect FFD signals to conBou";
  final parameter Integer kSurBou = kConBou + nConBou
    "Offset used to connect FFD signals to surBou";
  final parameter Integer kHeaPorAir = kSurBou + nSurBou
    "Offset used to connect FFD output signal to air heat port (to send average temperature from FFD to Modelica)";
//  final parameter Integer kUSha = kHeaPorAir + 1
  final parameter Integer kUSha = kSurBou + nSurBou
    "Offset used to connect FFD signals to input signal of shade";
  final parameter Integer kQRadAbs_flow = if haveShade then kUSha + nConExtWin else kUSha
    "Offset used to connect FFD signals to input signal that contains the radiation absorbed by the shade";
  // Because heaPorAir is only receiving T from FFD, but does not send Q_flow to FFD, there is no '+1' increment
  // for kTSha
  final parameter Integer kTSha = kHeaPorAir + 1
    "Offset used to connect FFD signals to output signal that contains the shade temperature";

  final parameter Integer kQConGai_flow = if haveShade then kQRadAbs_flow + nConExtWin else kQRadAbs_flow
    "Offset used to connect FFD signals to input signal for connect convective sensible heat gain";

  final parameter Integer kQLatGai_flow = kQConGai_flow + 1
    "Offset used to connect FFD signals to input signal for connect radiative heat gain";

  final parameter Integer kFluIntP = kQLatGai_flow + 1
    "Offset used to connect FFD signals to input signal for pressure from the fluid ports";

  final parameter Integer kFluIntM_flow = kFluIntP + 1
    "Offset used to connect FFD signals to input signals for mass flow rate from the fluid ports";
  final parameter Integer kFluIntT_inflow = kFluIntM_flow + nPorts
    "Offset used to connect FFD signals to input signals for inflowing temperature from the fluid ports";
  final parameter Integer kFluIntXi_inflow = kFluIntT_inflow + nPorts
    "Offset used to connect FFD signals to input signals for inflowing species concentration from the fluid ports";

  final parameter Integer kFluIntC_inflow = kFluIntXi_inflow + nPorts*Medium.nXi
    "Offset used to connect FFD signals to input signals for inflowing trace substances from the fluid ports";

  // Input signals to fluInt block
  final parameter Integer kFluIntT_outflow = if haveShade then kTSha + nConExtWin else kTSha
    "Offset used to connect FFD signals to outgoing temperature for the fluid ports";
  final parameter Integer kFluIntXi_outflow = kFluIntT_outflow+nPorts
    "Offset used to connect FFD signals to outgoing species concentration for the fluid ports";
  final parameter Integer kFluIntC_outflow = kFluIntXi_outflow+nPorts*Medium.nXi
    "Offset used to connect FFD signals to outgoing trace substances for the fluid ports";
  final parameter Integer kSen = kFluIntC_outflow+nPorts*Medium.nC
    "Offset used to connect FFD signals to output sensor";

  final parameter Integer nSur = kSurBou+nSurBou "Number of surfaces";
protected
  function assignSurfaceIdentifier

    input Integer nConExt(min=0) "Number of exterior constructions";
    input Integer nConExtWin(min=0) "Number of window constructions";
    input Integer nConPar(min=0) "Number of partition constructions";
    input Integer nConBou(min=0)
      "Number of constructions that have their outside surface exposed to the boundary of this room";
    input Integer nSurBou(min=0)
      "Number of surface heat transfer models that connect to constructions that are modeled outside of this room";
    input Integer nSur(min=1) "Total number of surfaces";

    input Boolean haveShade
      "Flag, set to true if any of the window in this room has a shade";
    /*
    // Declaration of counters used in the loop.
    // This could be computed (again) in this function, but using it
    // as a function arguments avoids code duplication.
    input Integer kConExt "Offset used to connect FFD signals to conExt";
    input Integer kConExtWin "Offset used to connect FFD signals to conExtWin";
    input Integer kGlaUns "Offset used to connect FFD signals to glaUns";
    input Integer kGlaSha "Offset used to connect FFD signals to glaSha";
    input Integer kConExtWinFra "Offset used to connect FFD signals to glaSha";
    input Integer kConPar_a "Offset used to connect FFD signals to conPar_a";
    input Integer kConPar_b "Offset used to connect FFD signals to conPar_b";
    input Integer kConBou "Offset used to connect FFD signals to conBou";
    input Integer kSurBou "Offset used to connect FFD signals to surBou";
*/
    // Declaration of construction data
    input String nameConExt[nConExt] "Surface name";
    input Modelica.SIunits.Area AConExt[nConExt] "Surface area";
    input Modelica.SIunits.Angle tilConExt[nConExt] "Surface tilt";
    input Buildings.Rooms.Types.CFDBoundaryConditions bouConConExt[nConExt]
      "Boundary condition";

    input String nameConExtWin[nConExtWin] "Surface name";
    input Modelica.SIunits.Area AConExtWin[nConExtWin] "Surface area";
    input Modelica.SIunits.Angle tilConExtWin[nConExtWin] "Surface tilt";
    input Buildings.Rooms.Types.CFDBoundaryConditions bouConConExtWin[nConExtWin]
      "Boundary condition";
    input Modelica.SIunits.Area AGla[nConExtWin] "Surface area";
    input Modelica.SIunits.Area AFra[nConExtWin] "Surface area";

    input String nameConPar[nConPar] "Surface name";
    input Modelica.SIunits.Area AConPar[nConPar] "Surface area";
    input Modelica.SIunits.Angle tilConPar[nConPar] "Surface tilt";
    input Buildings.Rooms.Types.CFDBoundaryConditions bouConConPar[nConPar]
      "Boundary condition";

    input String nameConBou[nConBou] "Surface name";
    input Modelica.SIunits.Area AConBou[nConBou] "Surface area";
    input Modelica.SIunits.Angle tilConBou[nConBou] "Surface tilt";
    input Buildings.Rooms.Types.CFDBoundaryConditions bouConConBou[nConBou]
      "Boundary condition";

    input String nameSurBou[nSurBou] "Surface name";
    input Modelica.SIunits.Area ASurBou[nSurBou] "Surface area";
    input Modelica.SIunits.Angle tilSurBou[nSurBou] "Surface tilt";
    input Buildings.Rooms.Types.CFDBoundaryConditions bouConSurBou[nSurBou]
      "Boundary condition";

    output FFDSurfaceIdentifier id[nSur] "Name of all surfaces";

  algorithm
      id := cat(1, {FFDSurfaceIdentifier(
                                  name=nameConExt[i],
                                  A=AConExt[i],
                                  til=tilConExt[i],
                                  bouCon=bouConConExt[i]) for i in 1:nConExt},
                   {FFDSurfaceIdentifier(
                                  name=nameConExtWin[i],
                                  A=AConExtWin[i],
                                  til=tilConExtWin[i],
                                  bouCon=bouConConExtWin[i]) for i in 1:nConExtWin},
                   {FFDSurfaceIdentifier(
                                  name=nameConExtWin[i] + " (glass, unshaded)",
                                  A=AGla[i],
                                  til=tilConExtWin[i],
                                  bouCon=bouConConExtWin[i]) for i in 1:nConExtWin},
                   {FFDSurfaceIdentifier(
                                  name=nameConExtWin[i] + " (glass, shaded)",
                                  A=AGla[i],
                                  til=tilConExtWin[i],
                                  bouCon=bouConConExtWin[i]) for i in 1:(if haveShade then nConExtWin else 0)},
                   {FFDSurfaceIdentifier(
                                  name=nameConExtWin[i] + " (frame)",
                                  A=AFra[i],
                                  til=tilConExtWin[i],
                                  bouCon=bouConConExtWin[i]) for i in 1:nConExtWin},
                   {FFDSurfaceIdentifier(
                                  name=nameConPar[i] + " (surface a)",
                                  A=AConPar[i],
                                  til=tilConPar[i],
                                  bouCon=bouConConPar[i]) for i in 1:nConPar},
                   {FFDSurfaceIdentifier(
                                  name=nameConPar[i] + " (surface b)",
                                  A=AConPar[i],
                                  til=tilConPar[i] + Modelica.Constants.pi/180,
                                  bouCon=bouConConPar[i]) for i in 1:nConPar},
                  {FFDSurfaceIdentifier(
                                  name=nameConBou[i],
                                  A=AConBou[i],
                                  til=tilConBou[i],
                                  bouCon=bouConConBou[i]) for i in 1:nConBou},
                  {FFDSurfaceIdentifier(
                                  name=nameSurBou[i],
                                  A=ASurBou[i],
                                  til=tilSurBou[i],
                                  bouCon=bouConSurBou[i]) for i in 1:nSurBou});
  end assignSurfaceIdentifier;

public
  Modelica.Blocks.Interfaces.RealInput QCon_flow
    "Convective sensible heat gains of the room"
    annotation (Placement(transformation(extent={{-280,-120},{-240,-80}})));
  Modelica.Blocks.Interfaces.RealInput QLat_flow
    "Latent heat gains for the room"
    annotation (Placement(transformation(extent={{-280,-180},{-240,-140}})));
  Modelica.Blocks.Math.Add QTotCon_flow
    "Total sensible convective heat flow rate added to the room"
    annotation (Placement(transformation(extent={{-160,-60},{-140,-40}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor senHeaFlo
    "Sensor for heat flow added through the port heaPorAir"
    annotation (Placement(transformation(extent={{-210,-10},{-190,10}})));
initial equation
   for i in 1:nPorts loop
     for j in 1:Medium.nXi loop
      Xi_outflow_fixed[(i-1)*Medium.nXi+j] = X_start[j];
     end for;
   end for;

  // Assignment of uStart
  for i in 1:kUSha loop
    uStart[i] = T_start;
  end for;
  if haveShade then
    for i in 1:nConExtWin loop
      uStart[kUSha+i] = 0;
      uStart[kQRadAbs_flow+i] = 0;
    end for;
  end if;
  uStart[kQConGai_flow+1] = 0;
  uStart[kQLatGai_flow+1] = 0;
  uStart[kFluIntP+1] = p_start;
  for i in 1:nPorts loop
    uStart[kFluIntM_flow+i] = 0;
    uStart[kFluIntT_inflow+i] = T_start;
    for j in 1:Medium.nXi loop
      uStart[kFluIntXi_inflow+(i-1)*Medium.nXi+j] =  X_start[j];
    end for;
    for j in 1:Medium.nC loop
      uStart[kFluIntC_inflow+(i-1)*Medium.nC+j] = C_start[j];
    end for;
  end for;

  // Assignment of yFixed
  for i in 1:kSurBou+nSurBou loop
    yFixed[i] = Q_flow_fixed[i];
  end for;
  yFixed[kHeaPorAir+1] = TRooAve_fixed;
  if haveShade then
    for i in 1:nConExtWin loop
      yFixed[kTSha+i] = TSha_fixed[i];
    end for;
  end if;
  for i in 1:nPorts loop
    yFixed[kFluIntT_outflow+i] = T_outflow_fixed[i];
    for j in 1:Medium.nXi loop
      yFixed[kFluIntXi_outflow+(i-1)*Medium.nXi+j] = Xi_outflow_fixed[(i-1)*Medium.nXi+j];
    end for;
    for j in 1:Medium.nC loop
      yFixed[kFluIntC_outflow+(i-1)*Medium.nC+j] = C_outflow_fixed[(i-1)*Medium.nC+j];
    end for;
  end for;
  for i in 1:nSen loop
    yFixed[kSen+i] = 0;
  end for;
equation
  //////////////////////////////////////////////////////////////////////
  // Data exchange with FFD block
  if haveConExt then
    for i in 1:nConExt loop
      if datConExt[i].boundaryCondition == Buildings.Rooms.Types.CFDBoundaryConditions.Temperature then
        connect(ffd.u[kConExt+i], ffdConExt[i].T_out) annotation (Line(
          points={{-42,190},{-60,190},{-60,224},{179,224}},
          color={0,0,127},
          smooth=Smooth.None));
        connect(ffd.y[kConExt+i], ffdConExt[i].Q_flow_in) annotation (Line(
          points={{-19,190},{60,190},{60,228},{178,228}},
          color={0,0,127},
          smooth=Smooth.None));
      else
        connect(ffd.u[kConExt+i], ffdConExt[i].T_in) annotation (Line(
          points={{-42,190},{-60,190},{-60,212},{179,212}},
          color={0,0,127},
          smooth=Smooth.None));
        connect(ffd.y[kConExt+i], ffdConExt[i].Q_flow_out) annotation (Line(
          points={{-19,190},{60,190},{60,216},{179,216}},
          color={0,0,127},
          smooth=Smooth.None));
      end if;
    end for;
  end if;

  if haveConExtWin then
    for i in 1:nConExtWin loop
      if datConExtWin[i].boundaryCondition == Buildings.Rooms.Types.CFDBoundaryConditions.Temperature then
        connect(ffd.u[kConExtWin+i], ffdConExtWin[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,172},{60,172},{60,184},{179,184}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConExtWin+i], ffdConExtWin[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,188},{178,188}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.u[kGlaUns+i], ffdGlaUns[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,124},{179,124}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kGlaUns+i], ffdGlaUns[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,128},{178,128}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.u[kConExtWinFra+i], ffdConExtWinFra[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,4},{179,4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConExtWinFra+i], ffdConExtWinFra[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,8},{178,8}},
            color={0,0,127},
            smooth=Smooth.None));
      else
        connect(ffd.u[kConExtWin+i], ffdConExtWin[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,172},{179,172}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConExtWin+i], ffdConExtWin[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,176},{179,176}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.u[kGlaUns+i], ffdGlaUns[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,112},{179,112}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kGlaUns+i], ffdGlaUns[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,116},{179,116}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.u[kConExtWinFra+i], ffdConExtWinFra[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-8},{179,-8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConExtWinFra+i], ffdConExtWinFra[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,-4},{179,-4}},
            color={0,0,127},
            smooth=Smooth.None));
      end if;
    end for;
  end if;

  if haveShade then
    for i in 1:nConExtWin loop
      if datConExtWin[i].boundaryCondition == Buildings.Rooms.Types.CFDBoundaryConditions.Temperature then
        connect(ffd.u[kGlaSha+i], ffdGlaSha[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,84},{179,84}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kGlaSha+i], ffdGlaSha[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,88},{178,88}},
            color={0,0,127},
            smooth=Smooth.None));
      else
        connect(ffd.u[kGlaSha+i], ffdGlaSha[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,72},{179,72}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kGlaSha+i], ffdGlaSha[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,76},{179,76}},
            color={0,0,127},
            smooth=Smooth.None));
      end if;
    end for;
  end if;

  if haveConPar then
    for i in 1:nConPar loop
      if datConPar[i].boundaryCondition == Buildings.Rooms.Types.CFDBoundaryConditions.Temperature then
        connect(ffd.u[kConPar_a+i], ffdConPar_a[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-56},{179,-56}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConPar_a+i], ffdConPar_a[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,-52},{178,-52}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.u[kConPar_b+i], ffdConPar_b[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-96},{179,-96}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConPar_b+i], ffdConPar_b[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,-92},{178,-92}},
            color={0,0,127},
            smooth=Smooth.None));
      else
        connect(ffd.u[kConPar_a+i], ffdConPar_a[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-68},{179,-68}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConPar_a+i], ffdConPar_a[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,-64},{179,-64}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.u[kConPar_b+i], ffdConPar_b[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-108},{179,-108}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConPar_b+i], ffdConPar_b[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,-104},{179,-104}},
            color={0,0,127},
            smooth=Smooth.None));
      end if;
    end for;
  end if;

  if haveConBou then
    for i in 1:nConBou loop
      if datConBou[i].boundaryCondition == Buildings.Rooms.Types.CFDBoundaryConditions.Temperature then
        connect(ffd.u[kConBou+i], ffdConBou[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-156},{179,-156}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConBou+i], ffdConBou[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,-152},{178,-152}},
            color={0,0,127},
            smooth=Smooth.None));
      else
        connect(ffd.u[kConBou+i], ffdConBou[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-168},{179,-168}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kConBou+i], ffdConBou[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,-164},{179,-164}},
            color={0,0,127},
            smooth=Smooth.None));
      end if;
    end for;
  end if;

  if haveSurBou then
    for i in 1:nSurBou loop
      if surBou[i].boundaryCondition == Buildings.Rooms.Types.CFDBoundaryConditions.Temperature then
        connect(ffd.u[kSurBou+i], ffdSurBou[i].T_out)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-216},{179,-216}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kSurBou+i], ffdSurBou[i].Q_flow_in)
            annotation (Line(
            points={{-19,190},{60,190},{60,-212},{178,-212}},
            color={0,0,127},
            smooth=Smooth.None));
      else
        connect(ffd.u[kSurBou+i], ffdSurBou[i].T_in)
            annotation (Line(
            points={{-42,190},{-60,190},{-60,-228},{179,-228}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(ffd.y[kSurBou+i], ffdSurBou[i].Q_flow_out)
            annotation (Line(
            points={{-19,190},{60,190},{60,-224},{179,-224}},
            color={0,0,127},
            smooth=Smooth.None));
      end if;
    end for;
  end if;

  connect(ffdConExt.port, conExt) annotation (Line(
      points={{200,220},{240,220}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdConExtWin.port, conExtWin) annotation (Line(
      points={{200,180},{240,180}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdGlaUns.port, glaUns) annotation (Line(
      points={{200,120},{220,120},{220,120},{240,120}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdGlaSha.port, glaSha) annotation (Line(
      points={{200,80},{240,80}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdConExtWinFra.port, conExtWinFra) annotation (Line(
      points={{200,0},{242,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdConPar_a.port, conPar_a) annotation (Line(
      points={{200,-60},{242,-60}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdConPar_b.port, conPar_b) annotation (Line(
      points={{200,-100},{242,-100}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdConBou.port, conBou) annotation (Line(
      points={{200,-160},{222,-160},{222,-160},{242,-160}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ffdSurBou.port, conSurBou) annotation (Line(
      points={{200,-220},{241,-220}},
      color={191,0,0},
      smooth=Smooth.None));
  // Connections to heat port of air volume
  connect(ffd.y[kHeaPorAir+1], ffdHeaPorAir.T) annotation (Line(
      points={{-19,190},{60,190},{60,8.88178e-16},{-138,8.88178e-16}},
      color={0,0,127},
      smooth=Smooth.None));
  // Connections to shade
  if haveShade then
    for i in 1:nConExtWin loop
      connect(ffd.u[kUSha+i], uSha[i]) annotation (Line(
          points={{-42,190},{-60,190},{-60,200},{-260,200}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ffd.u[kQRadAbs_flow+i], QRadAbs_flow[i]) annotation (Line(
          points={{-42,190},{-60,190},{-60,90},{-260,90}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ffd.y[kTSha+i], TSha[i]) annotation (Line(
          points={{-19,190},{60,190},{60,60},{-250,60}},
          color={0,0,127},
          smooth=Smooth.None));
    end for;
  end if;

  // Connection for heat gain that is added to the room
  // (averaged over the whole room air volume)
  connect(QTotCon_flow.y, ffd.u[kQConGai_flow+1]) annotation (Line(
      points={{-139,-50},{-60,-50},{-60,190},{-42,190}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(QLat_flow, ffd.u[kQLatGai_flow+1]) annotation (Line(
      points={{-260,-160},{-190,-160},{-190,-72},{-60,-72},{-60,190},{-42,190},{
          -42,190}},
      color={0,0,127},
      smooth=Smooth.None));

  // Connections to fluid port
  connect(ports, fluInt.ports) annotation (Line(
      points={{0,-238},{0,-198}},
      color={0,127,255},
      smooth=Smooth.None));

  // Output signals from fluInt block

  // The pressure of the air volume will be sent from Modelica to FFD
  connect(ffd.u[kFluIntP+1], fluInt.p) annotation (Line(
      points={{-42,190},{-60,190},{-60,-180},{-11,-180}},
      color={0,0,127},
      smooth=Smooth.None));
  for i in 1:nPorts loop
    connect(ffd.u[kFluIntM_flow+i], fluInt.m_flow[i]) annotation (Line(
        points={{-42,190},{-60,190},{-60,-184},{-11,-184}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(ffd.u[kFluIntT_inflow+i], fluInt.T_inflow[i]) annotation (Line(
        points={{-42,190},{-60,190},{-60,-188},{-11,-188}},
        color={0,0,127},
        smooth=Smooth.None));
    for j in 1:Medium.nXi loop
      connect(ffd.u[kFluIntXi_inflow+(i-1)*Medium.nXi+j], fluInt.Xi_inflow[(i-1)*Medium.nXi+j]) annotation (Line(
          points={{-42,190},{-60,190},{-60,-192},{-11,-192}},
          color={0,0,127},
          smooth=Smooth.None));
    end for;
    for j in 1:Medium.nC loop
      connect(ffd.u[kFluIntC_inflow+(i-1)*Medium.nC+j], fluInt.C_inflow[(i-1)*Medium.nC+j]) annotation (Line(
          points={{-42,190},{-60,190},{-60,-196},{-11,-196}},
          color={0,0,127},
          smooth=Smooth.None));
    end for;
  end for;
  // Input signals to fluInt block
  // The pressures of ports[2:nPorts] will be sent from FFD to Modelica
  for i in 1:nPorts loop
    connect(ffd.y[kFluIntT_outflow+i], fluInt.T_outflow[i]) annotation (Line(
        points={{-19,190},{60,190},{60,-188},{12,-188}},
        color={0,0,127},
        smooth=Smooth.None));
    for j in 1:Medium.nXi loop
      connect(ffd.y[kFluIntXi_outflow+(i-1)*Medium.nXi+j], fluInt.Xi_outflow[(i-1)*Medium.nXi+j]) annotation (Line(
          points={{-19,190},{60,190},{60,-192},{12,-192}},
          color={0,0,127},
          smooth=Smooth.None));
    end for;
    for j in 1:Medium.nC loop
      connect(ffd.y[kFluIntC_outflow+(i-1)*Medium.nC+j], fluInt.C_outflow[(i-1)*Medium.nC+j]) annotation (Line(
          points={{-19,190},{60,190},{60,-196},{12,-196}},
          color={0,0,127},
          smooth=Smooth.None));
    end for;
  end for;

  // Connections for sensor signal
  if haveSensor then
    for i in 1:nSen loop
      connect(ffd.y[kSen+i], yCFD[i]) annotation (Line(
        points={{-19,190},{60,190},{60,-234},{180,-234},{180,-250}},
        color={0,0,127},
        smooth=Smooth.None));
    end for;
  end if;

  connect(heaPorAir, senHeaFlo.port_a) annotation (Line(
      points={{-240,0},{-210,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(senHeaFlo.port_b, ffdHeaPorAir.port) annotation (Line(
      points={{-190,0},{-160,0}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(senHeaFlo.Q_flow, QTotCon_flow.u1) annotation (Line(
      points={{-200,-10},{-200,-46},{-162,-46},{-162,-44}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(QCon_flow, QTotCon_flow.u2) annotation (Line(
      points={{-260,-100},{-200,-100},{-200,-56},{-162,-56}},
      color={0,0,127},
      smooth=Smooth.None));

  annotation (
    preferredView="info",
    Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-240,-240},{240,
            240}}), graphics),
    Icon(coordinateSystem(preserveAspectRatio=false,extent={{-240,-240},{240,240}}),
                    graphics={
          Rectangle(
          extent={{-144,184},{148,-200}},
          pattern=LinePattern.None,
          lineColor={0,0,0},
          fillColor={170,213,255},
          fillPattern=FillPattern.Sphere)}),
    Documentation(info="<html>
<p>
This model computes the heat and mass balance of the air using Fast Fluid Dynamics.
</p>
<h4>Conventions</h4>
<p>
The following conventions are made:
</p>
<ol>
<li>
<p>
The port <code>heaPorAir</code> contains the average room air temperature, defined as
</p>
<p align=\"center\" style=\"font-style:italic;\">
  T<sub>a</sub> = 1 &frasl; V &nbsp; &int;<sub>V</sub> T(dV) &nbsp; dV,
</p>
<p>
where <i>T<sub>a</sub></i> is the average room air temperature, <i>V</i> is the room air volume
and <i>T(dV)</i> is the room air temperature in the control volume <i>dV</i>.
The average room air temperature <i>T<sub>a</sub></i> is computed by the FFD program.
</p>
</li>
<li>
If a model injects heat to <code>heaPorAir</code>, then the heat will be distributed to all
cells. The amount of heat flow rate that each cell exchanges with <code>heaPorAir</code> is
proportional to its volume.
</li>
<li>
If a construction is not present, or if no shade is present, or if no air stream is connected to <code>ports</code>,
then no variables are exchanged for this quantity with the block <code>ffd</code>.
</li>
<li>
The variables for the heat ports that connect to wall surfaces
are the heat flow rate and the temperature.
The interface with the block <code>ffd</code> is done through the instances
of the model
<a href=\"modelica://Buildings.Rooms.BaseClasses.FFDSurfaceInterface\">
Buildings.Rooms.BaseClasses.FFDSurfaceInterface</a>.
This model has four ports, and depepending on the parameter
<code>bouCon</code>, two of these ports are conditionally removed.
This allows to use the parameter <code>bouCon</code> to specify whether
the surface should be used with a temperature or a heat flow rate
boundary condition.
Therefore, the inputs and outputs to the instance <code>ffd</code>
are either temperatures or heat flow rates.
The parameter <code>surIde</code> of this model, which is also
propagated to the instance <code>ffd</code>, declares what
type of boundary condition is used.
</li>
<li>
The variables of the connector <code>ports</code> are exchanged with the FFD block
through the instance <code>intFlu</code>. Its output and input signals are as follows:
<ul>
<li>
Input to the FFD block is a vector <code>[p, m_flow[nPorts], T_inflow[nPorts], 
X_inflow[nPorts*Medium.nXi], C_inflow[nPorts*Medium.nC]]</code>.
The quantity <code>p</code> is the total pressure of the fluid ports (all fluid ports have the same
total pressure). 
Therefore, the flow resistance of the diffusor or exhaust grill must be computed in the 
Modelica HVAC system that is connected to the room model.
The quantities <code>X_inflow</code> and <code>C_inflow</code> (or <code>X_inflow</code> and <code>C_inflow</code>)
are vectors with components <code>X_inflow[1:Medium.nXi]</code> and <code>C_inflow[1:Medium.nC]</code>.
For example, for moist air, <code>X_inflow</code> has one element which is equal to the mass fraction of air,
relative to the total air mass and not the dry air.
</li>
<li>
Output from the FFD block is a vector 
<code>[T_outflow[nPorts], X_outflow[nPorts*Medium.nXi], C_outflow[nPorts*Medium.nC]]</code>.
The quantities <code>*_outflow</code> are the fluid properties of the cell to which the port is
connected. 
</li>
<li>
If <code>Medium.nXi=0</code> (e.g., for dry air) or <code>Medium.nC=0</code>, then these signals are not present as input/output signals of the FFD block.
</li>
</ul>
The quantities that are exchanged between the programs are defined as follows:
<ul>
<li>
For the mass flow rate of the fluid port, 
we exchange <i>m<sub>e</sub> = 1 &frasl; &Delta; t &int;<sub>&Delta; t</sub> m(s) dt</i>.
</li>
<li>
For the temperature, species concentration and trace substances of the fluid port, we exchange 
<i>X = 1 &frasl; (m<sub>e</sub> &nbsp; &Delta; t) &int;<sub>&Delta; t</sub> m(s) &nbsp; X(s) dt</i>.
Note that for the first implementation, FFD does only compute a bulk mass balance for <code>Xi</code>.
It does not do a moisture balance for each cell.
However, for trace substances <code>C</code>, FFD does a contaminant balance for each cell
and return <code>C_outflow</code> to be the contaminant concentration of that cell.
</li>
<li>
For the surface temperatures, 
we exchange <i>T<sub>e</sub> = 1 &frasl; &Delta; t &int;<sub>&Delta; t</sub> T(s) dt</i>.
</li>
<li>
For the surface heat flow rates, 
we exchange <i>Q<sub>e</sub> = 1 &frasl; &Delta; t &int;<sub>&Delta; t</sub> Q(s) dt</i>.
</li>
</ul>
</li>
</ol>
<p>
The data exchange with the CFD interface is done through the instance
<code>ffd</code>, and implemented in
<a href=\"modelica://Buildings.Rooms.BaseClasses.FFDExchange\">
Buildings.Rooms.BaseClasses.FFDExchange</a>.
This block exchanges the following data with the CFD simulation:
</p>
<ul>
<li>
During the initialzation, the following data are sent from Modelica to CFD:
<ul>
<li>
An array of strings where each element is the name of the surface, 
as declared by
the user when instantiating the model
<a href=\"Buildings.Rooms.FFD\">Buildings.Rooms.FFD</a>.
Let us call this array <code>name</code>.
The orders of elements in this array are as follows:
<ol>
<li>
The first 
<code>nConExt</code> elements are the names of the exterior constructions 
declared as <code>datConExt</code>. The order is the same as 
in the declaration of <code>datConExt</code>.
</li>
<li>
<code>nConExtWin</code> elements are the names of the exterior constructions 
declared as <code>datConExtWin</code>. These constructions embed windows
and a frame. Therefore, what follows are
<code>nConExtWin</code> elements where each string is the same as above, 
but <code>' (glass, unshaded)'</code> has been appended, 
then -- if and only if the window has a shade --
<code>nConExtWin</code> elements follow with 
<code>' (glass, shaded)'</code> appended, and,
finally,
<code>nConExtWin</code> elements follow with <code>' (frame)'</code> appended.
</li>
<li>
<code>nConPar</code> elements for the surface <code>a</code> of <code>datConPar</code>.
To these names, the string <code>' (surface a)'</code> is appended.
Next, there are <code>nConPar</code> elements with  <code>' (surface b)'</code> appended.
</li>
<li>
<code>nConBou</code> elements for the surfaces of <code>datConBou</code>.
</li>
<li>
<code>nSurBou</code> elements for the surfaces of <code>nSurBou</code>.
</li>
</ol>
</li>
<li>
Using the same order, there is also an array for the areas of the surfaces <code>A</code>,
an array for the surface tilt <code>til</code>
and the type of the boundary conditions <code>bouCon</code> for each of these surfaces.
If <code>bouCon[i] = 1</code>, 
then temperature is sent from Modelica to CFD.
If <code>bouCon[i] = 2</code>, then
heat flow rate is sent from Modelica to CFD.
</li>
</ul>
</li>
<li>
During the time integration, and array <code>u</code> is sent from Modelica to CFD, and Modelica
receives an array <code>y</code> from CFD.
The elements of the array <code>u</code> are as follows:
<ol>
<li>
Either temperature or heat flow rate boundary conditions, 
in the same order as the array <code>name</code>. The units are <i>[K]</i> or <i>[W]</i>.
The array <code>bouCon</code> that is sent during the 
initialization declares the type of boundary
condition.
There are <code>nSur</code> elements for surfaces.
</li>
<li>
If at least one window in the room has a shade, then the next 
<code>nConExtWin</code>
elements are the shading control signals. <code>u=0</code> means 
that the shade is not deployed,
and <code>u=1</code> means that the shade is 
completely deployed (blocking solar radiation).
If there is no window in the room, then these elements are not present.
</li>
<li>
If at least one window in the room has a shade, then the next <code>nConExtWin</code>
elements are the radiations in <i>[W]</i> that are absorbed by the 
respective shades.
If there is no window in the room, then these elements are not present.
</li>
<li>
The convective sensible heat input into the room in <i>[W]</i>, which is a scalar.
A positive value means that heat is added to the room.
</li>
<li>
The latent heat input into the room in <i>[W]</i>, which is a scalar.
A positive value means that moisture is added to the room.
</li>
<li>
The next element is the room average static pressure in <i>[Pa]</i>.
</li>
<li>
The next <code>nPorts</code> elements are the mass flow rates
into the room in <i>[kg/s]</i>. 
A positive value is used if the air flows into the room,
otherwise the value is negative.
The first element is connected to <code>ports[1]</code>, the second to
<code>ports[2]</code> etc.
</li>
<li>
The next <code>nPorts</code> elements are the air temperatures 
that the medium has
<i>if it were flowing into the room</i>, e.g., the \"inflowing medium\" 
computed based on <code>inStream(h_outflow)</code>.
</li>
<li>
The next <code>nPorts*Medium.nXi</code> elements are the 
species concentration of the inflowing
medium. 
The first <code>Medium.nXi</code> elements are for port <i>1</i>, then for 
port <i>2</i> etc.
The units are in <i>[kg/kg]</i> total mass, and not in  <i>[kg/kg]</i> dry air.
</li>
<li>
The next <code>nPorts*Medium.nC</code> elements are the trace substances 
of the inflowing
medium. 
The first <code>Medium.nC</code> elements are for port <i>1</i>, then for 
port <i>2</i> etc.
</li>
</ol>
The elements of the array <code>y</code> that is sent from CFD to Modelica are as follows:
<ol>
<li>
Either temperature or heat flow rate at the surfaces,
in the same order as the array <code>name</code>.
The array <code>bouCon</code> that is sent during the 
initialization declares the type of boundary
condition.
If <code>bouCon[i] = 1</code>, then heat flow rate in <i>[W]</i>
is sent from CFD to Modelica.
If <code>bouCon[i] = 2</code>, then temperature in <i>[K]</i>
is sent from CFD to Modelica.
There are <code>nSur</code> elements for surfaces.
</li>
<li>
The average room air temperature in <i>[K]</i>.
</li>
<li>
If the room has at least one window with a shade, then the next
<code>nConExtWin</code> elements are the temperature of the
shade in <i>[K]</i>.
</li>
<li>
The next <code>nPorts</code> elements are the air temperatures
in <i>[K]</i>
of the cells that are connected to the inlet or outlet diffusor 
of <code>ports[1], ports[2], etc.</code>.
</li>
<li>
The next <code>nPorts*Medium.nXi</code> elements are the 
species concentration of the cells to which the ports are 
connected.
The first <code>Medium.nXi</code> elements are for port <i>1</i>, then for 
port <i>2</i> etc.
The units are in <i>[kg/kg]</i> total mass, and not in <i>[kg/kg]</i> dry air.
</li>
<li>
The next <code>nPorts*Medium.nC</code> elements are the trace substances 
of the cells to which the ports are connected to.
The first <code>Medium.nC</code> elements are for port <i>1</i>, then for 
port <i>2</i> etc.
</li>
</ol>
</li>
</ul>
</html>",
revisions="<html>
<ul>
<li>
July 17, 2013, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
end AirHeatMassBalanceFFD;
