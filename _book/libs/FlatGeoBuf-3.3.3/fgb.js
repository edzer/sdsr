LeafletWidget.methods.addFlatGeoBuf = function (layerId,
                                                group,
                                                url,
                                                popup,
                                                label,
                                                style,
                                                options,
                                                className,
                                                scale,
                                                scaleFields) {

  var map = this;
  var gl = false;
  var pane;

  if (options.pane === undefined) {
    pane = 'overlayPane';
  } else {
    pane = options.pane;
  }

  var data_fl = document.getElementById(layerId + '-1-attachment');

  if (data_fl === null) {
    data_fl = url;
  } else {
    data_fl = data_fl.href;
  }

  var popUp;
  var colnames = [];

  function handleHeaderMeta(headerMeta) {
    headerMeta.columns.forEach(function(col) {
      colnames.push(col.name);
    });
  }

  function handleResponse(response) {
    // use fgb JavaScript API to iterate stream into results (features as geojson)
    // NOTE: would be more efficient with a special purpose Leaflet deserializer
    let it = flatgeobuf.deserialize(response.body, undefined, handleHeaderMeta);
    // handle result
    function handleResult(result) {
        if (!result.done) {
          if (gl) {
            map.layerManager.addLayer(
              L.glify.shapes({
                map: map,
                data: result.value,
                className: group
              }).glLayer, null, null, group);
            it.next().then(handleResult);
          } else {

            if (popup) {
              pop = makePopup(popup, className);
            } else {
              pop = null;
            }

            if (scaleFields === undefined &
                result.value.properties !== undefined) {
              var vls = Object.values(style);
              scaleFields = [];
              vls.forEach(function(name) {
                //if (name in colnames) {
                if (colnames.includes(name)) {
                  scaleFields.push(true);
                } else {
                  scaleFields.push(false);
                }
              });
            }

            lyr = L.geoJSON(result.value, {
              pointToLayer: function (feature, latlng) {
                  return L.circleMarker(latlng, options);
              },
              style: function(feature) {
                return updateStyle(style, feature, scale, scaleFields);
              },
              onEachFeature: pop,
              pane: pane
            });

            if (label) {
              if (Object.keys(result.value.properties).includes(label)) {
                lyr.bindTooltip(function (layer) {
                   return layer.feature.properties[label].toString();
                }, {sticky: true});
              }
            }

            map.layerManager.addLayer(lyr, null, null, group);
            it.next().then(handleResult);
          }
        }
    }
    it.next().then(handleResult);
  }

  fetch(data_fl) //, {mode: 'no-cors'})
  .then(handleResponse);

  //map.fitBounds(lyr.getBounds());
  //map.layerManager.addLayer(layer, null, null, group);
};

function makePopup(popup, className) {
  if (popup === true) {
    pop = function(feature, layer) {
      popUp = json2table(feature.properties, className);
      layer.bindPopup(popUp, { maxWidth: 2000 });
    };
  } else if (typeof(popup) === "string") {
    pop = function(feature, layer) {
      if (feature.properties !== undefined && popup in feature.properties) {
        popup = popup.split();
        popUp = json2table(
          pick(feature.properties, popup),
          className
        );
      } else {
        popUp = popup;
      }
      layer.bindPopup(popUp, { maxWidth: 2000 });
    };
  } else {
    pop = function(feature, layer) {
      popUp = json2table(
        pick(feature.properties, popup),
        className
      );
      layer.bindPopup(popUp, { maxWidth: 2000 });
    };
  }
  return pop;
}


function json2table(json, cls) {
  var cols = Object.keys(json);
  var vals = Object.values(json);

  var tab = "";

  for (let i = 0; i < cols.length; i++) {
    tab += "<tr><th>" + cols[i] + "&emsp;</th>" +
    "<td align='right'>" + vals[i] + "&emsp;</td></tr>";
  }

  return "<table class=" + cls + ">" + tab + "</table>";

}


/**
 * from https://gomakethings.com/how-to-create-a-new-object-with-only-a-subject-of-properties-using-vanilla-js/
 *
 *
 * Create a new object composed of properties picked from another object
 * (c) 2018 Chris Ferdinandi, MIT License, https://gomakethings.com
 * @param  {Object} obj   The object to pick properties from
 * @param  {Array}  props An array of properties to use
 * @return {Object}       The new object
 */
function pick(obj, props) {

	'use strict';

	// Make sure object and properties are provided
	if (!obj || !props) return;

	// Create new object
	var picked = {};

	// Loop through props and push to new object
	props.forEach(function(prop) {
		picked[prop] = obj[prop];
	});

	// Return new object
	return picked;

}


function updateStyle(style_obj, feature, scale, scaleValues) {
  var cols = Object.keys(style_obj);
  var vals = Object.values(style_obj);

  var out = {};

  for (let i = 0; i < cols.length; i++) {
    if (vals[i] === null) {
      out[cols[i]] = feature.properties[cols[i]];
    } else {
      if (scaleValues !== undefined) {
        //if (Object.keys(feature.properties).includes(vals[i])) {
        if (scaleValues[i] === true) {
          vals[i] = rescale(
            feature.properties[vals[i]]
            , scale[cols[i]].to[0]
            , scale[cols[i]].to[1]
            , scale[cols[i]].from[0]
            , scale[cols[i]].from[1]
          );
        }
      }
      out[cols[i]] = vals[i];
    }
  }

  return out;
}


function rescale(value, to_min, to_max, from_min, from_max) {
  if (value === undefined) {
    value = from_min;
  }
  return (value - from_min) / (from_max - from_min) * (to_max - to_min) + to_min;
}
