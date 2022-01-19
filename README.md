# GeoTransSwift
Swift 좌표계 변환 코드.

Java 좌표계 변환 코드 그대로 포팅했습니다.

Java 좌표계 변환 코드 링크: https://www.androidpub.com/android_dev_info/1043970

# 샘플코드:
        let in_pt = GeoPoint(longitude: 127.024612, latitude: 37.532600)
        
        let tm_pt = GeoTrans.shared.convert(srctype: .wgs84, dsttype: .tm, in_pt: in_pt)
        print("tm_pt \(tm_pt.longitude) \(tm_pt.latitude)")
        
        let katec_pt = GeoTrans.shared.convert(srctype: .tm, dsttype: .katec, in_pt: tm_pt)
        print("katec_pt \(katec_pt.longitude) \(katec_pt.latitude)")
        
        let out_pt = GeoTrans.shared.convert(srctype: .katec, dsttype: .wgs84, in_pt: katec_pt)
        print("out_pt \(out_pt.longitude) \(out_pt.latitude)")
        
        let in2_pt = GeoPoint(longitude: 128.0, latitude: 38.0)
        let distance = GeoTrans.shared.getDistancebyGeo(in_pt, in2_pt)
        print("distance \(distance)")
        
        let grsPoint = GeoTrans.shared.convert(srctype: .wgs84, dsttype: .grs80, in_pt: in_pt)
        print("grsPoint \(grsPoint.longitude) \(grsPoint.latitude)")
