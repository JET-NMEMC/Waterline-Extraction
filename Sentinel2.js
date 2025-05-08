var S2TOA = ee.ImageCollection("COPERNICUS/S2_HARMONIZED");
var S2SR = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');

var scaleImage = function (image, S2type) {
    var opticalBands = (S2type === 'S2TOA')
        ? image.select('B.*').divide(10000)
        : image.select('B.*');
    return image.addBands(opticalBands, null, true);
};

// Reference
// Normalized Burn Ratio, Elvidge, production
// Hammer, Kraft, and Steele (Data Lab at WRI)
// GFW-Fires, prototype
// https://gist.github.com/robinkraft/077c14d35a50a8b31581
var sentinelCloudScore = function (img, roi) {

    var rescale = function (img, thresholds) {
        return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
    };

    var score = ee.Image(1);
    // Clouds are reasonably bright in the blue band.
    score = score.min(rescale(img.select(['B2']), [0.1, 0.3]));
    // Clouds are reasonably bright in all visible bands.
    score = score.min(rescale(img.select(['B4']).add(img.select(['B3'])).add(img.select('B2')), [0.2, 0.8]));
    // Clouds are reasonably bright in all infrared bands.
    score = score.min(rescale(img.select(['B8']).add(img.select(['B11'])).add(img.select('B12')), [0.3, 0.8]));
    // However, clouds are not snow.
    var ndsi = img.normalizedDifference(['B3', 'B11']);
    score = score.min(rescale(ndsi, [0.8, 0.6])).rename(['cloudScore']);

    var iscloud = score.gt(0.4).multiply(100).rename(['cloud']);

    var cloudiness = iscloud.reduceRegion({
        reducer: 'mean',
        geometry: roi,
        scale: 30,
    });

    return img.addBands(score).addBands(iscloud).set(cloudiness);
};



var addIndices = function (image) {
    return image.addBands(image.normalizedDifference(['B3', 'B8']).float().rename('ndwi'))
        .addBands(image.normalizedDifference(['B3', 'B11']).float().rename('mndwi'))
        .addBands(image.normalizedDifference(['B8', 'B4']).float().rename('ndvi'));
};



var S2 = function (S2type, startdate, enddate, startmonth, endmonth, point, roi, cloudmax1, cloudmax2) {
    var datainput = (S2type === 'S2TOA') ? S2TOA : S2SR;
    var col = datainput.filterDate(startdate, enddate)
        .filter(ee.Filter.calendarRange(startmonth, endmonth, "month"))
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudmax1))
        .filterBounds(point)
        .filter(ee.Filter.eq('PRODUCTION_DEM_TYPE', null))
        .map(function (image) { return scaleImage(image, S2type); })
        .map(addIndices)
        .sort("system:time_start");

    if (cloudmax2 > 0) {
        col = col.map(function (img) { return sentinelCloudScore(img, roi); });
        var c1 = col.filter(ee.Filter.lt('cloud', cloudmax2));
        var c2 = col.filter(ee.Filter.gt('cloud', cloudmax2));
        return { 'c1': c1, 'c2': c2 };
    } else {
        return { 'c1': col, 'c2': null };
    }
};

exports.S2TOA = function (startdate, enddate, startmonth, endmonth, point, roi, cloudmax1, cloudmax2) {
    return S2('S2TOA', startdate, enddate, startmonth, endmonth, point, roi, cloudmax1, cloudmax2);
};

exports.S2SR = function (startdate, enddate, startmonth, endmonth, point, roi, cloudmax1, cloudmax2) {
    return S2('S2SR', startdate, enddate, startmonth, endmonth, point, roi, cloudmax1, cloudmax2);
};
