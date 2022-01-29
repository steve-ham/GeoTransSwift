# GeoTransSwift
Swift 좌표계 변환 코드.

Java 좌표계 변환 코드 그대로 포팅했습니다.

Java 좌표계 변환 코드 링크: https://www.androidpub.com/android_dev_info/1043970

# 샘플코드:
        let wgs84GeoPoint = GeoPoint(longitude: 127.024612, latitude: 37.532600)
        
        let tmGeoPoint = GeoTrans.shared.convert(from: .wgs84, to: .tm, geoPoint: wgs84GeoPoint)
        print("tmGeoPoint \(tmGeoPoint.longitude) \(tmGeoPoint.latitude)")
        
        let katecGeoPoint = GeoTrans.shared.convert(from: .tm, to: .katec, geoPoint: tmGeoPoint)
        print("katecGeoPoint \(katecGeoPoint.longitude) \(katecGeoPoint.latitude)")
        
        let wgs84GeoPoint2 = GeoTrans.shared.convert(from: .katec, to: .wgs84, geoPoint: katecGeoPoint)
        print("wgs84GeoPoint2 \(wgs84GeoPoint2.longitude) \(wgs84GeoPoint2.latitude)")
        
        let wgs84GeoPoint3 = GeoPoint(longitude: 128.0, latitude: 38.0)
        let distance = GeoTrans.shared.getDistancebyGeo(wgs84GeoPoint, wgs84GeoPoint3)
        print("distance \(distance)")
        
        let grs80GeoPoint = GeoTrans.shared.convert(from: .wgs84, to: .grs80, geoPoint: wgs84GeoPoint)
        print("grs80GeoPoint \(grs80GeoPoint.longitude) \(grs80GeoPoint.latitude)")
