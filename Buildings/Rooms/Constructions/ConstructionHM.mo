within Buildings.Rooms.Constructions;
model ConstructionHM "Model for an opaque construction that has no window"
  extends Buildings.Rooms.Constructions.BaseClasses.PartialConstructionHM(
    final AOpa=A);

  annotation (
defaultComponentName="conOpa",
Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-300,-300},
            {300,300}},
        initialScale=0.1)),
                          Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-300,-300},{300,300}},
        initialScale=0.1), graphics={
        Rectangle(
          extent={{2,260},{62,60}},
          lineColor={0,0,0},
          fillColor={97,97,135},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-54,260},{2,60}},
          lineColor={0,0,0},
          fillColor={183,123,103},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-78,260},{-54,60}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-54,260},{-68,274},{-16,274},{2,260},{-54,260}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={183,123,103},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-78,260},{-86,274},{-68,274},{-54,260},{-78,260}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{2,260},{-16,274},{42,274},{62,260},{2,260}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={97,97,135},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-86,274},{-86,82},{-78,60},{-78,260},{-86,274}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid)}),
    Documentation(
    info="<html>
    This model is used to compute heat and moisture transfer through opaque constructions inside the 
room model.
The model uses the record <code>layers</code> to access the material properties
of the opaque construction. The heat and moisture transfers are computed in the instance
<code>opa</code>, which uses the model 
<a href=\"modelica://Buildings.HeatTransfer.Conduction.MultiLayerHM\">
Buildings.HeatTransfer.Conduction.MultiLayerHM</a>.
</html>",
revisions="<html>
<ul>
<li>
July 6 2014, by Antoine Hudzik:<br/>
First implementation.
</li>
</ul>
</html>"));
end ConstructionHM;
