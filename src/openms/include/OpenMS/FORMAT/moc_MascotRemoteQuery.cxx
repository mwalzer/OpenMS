/****************************************************************************
** Meta object code from reading C++ file 'MascotRemoteQuery.h'
**
** Created: Mon Mar 3 21:41:29 2014
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MascotRemoteQuery.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MascotRemoteQuery.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_OpenMS__MascotRemoteQuery[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      41,   40,   40,   26, 0x05,
      53,   48,   40,   26, 0x05,

 // slots: signature, parameters, type, tag, flags
      86,   40,   40,   26, 0x0a,
      92,   40,   40,   26, 0x08,
     103,   48,   40,   26, 0x08,
     155,  138,   40,   26, 0x08,
     208,  185,   40,   26, 0x08,
     261,  238,   40,   26, 0x08,
     302,  291,   40,   26, 0x08,
     332,  326,   40,   26, 0x08,
     360,  354,   40,   26, 0x08,
     391,  375,   40,   26, 0x08,
     431,   40,   40,   26, 0x08,
     439,   40,   40,   26, 0x08,
     464,  451,   40,   26, 0x08,
     484,   48,   40,   26, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_OpenMS__MascotRemoteQuery[] = {
    "OpenMS::MascotRemoteQuery\0OPENMS_DLLAPI\0"
    "\0done()\0resp\0gotRedirect(QHttpResponseHeader)\0"
    "run()\0timedOut()\0readyReadSlot(QHttpResponseHeader)\0"
    "request_id,error\0httpRequestFinished(int,bool)\0"
    "bytes_read,bytes_total\0"
    "httpDataReadProgress(int,int)\0"
    "bytes_sent,bytes_total\0"
    "httpDataSendProgress(int,int)\0request_id\0"
    "httpRequestStarted(int)\0state\0"
    "httpStateChanged(int)\0error\0httpDone(bool)\0"
    "response_header\0readResponseHeader(QHttpResponseHeader)\0"
    "login()\0execQuery()\0results_path\0"
    "getResults(QString)\0"
    "followRedirect(QHttpResponseHeader)\0"
};

const QMetaObject OpenMS::MascotRemoteQuery::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_OpenMS__MascotRemoteQuery,
      qt_meta_data_OpenMS__MascotRemoteQuery, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &OpenMS::MascotRemoteQuery::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *OpenMS::MascotRemoteQuery::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *OpenMS::MascotRemoteQuery::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_OpenMS__MascotRemoteQuery))
        return static_cast<void*>(const_cast< MascotRemoteQuery*>(this));
    if (!strcmp(_clname, "DefaultParamHandler"))
        return static_cast< DefaultParamHandler*>(const_cast< MascotRemoteQuery*>(this));
    return QObject::qt_metacast(_clname);
}

int OpenMS::MascotRemoteQuery::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: done(); break;
        case 1: gotRedirect((*reinterpret_cast< const QHttpResponseHeader(*)>(_a[1]))); break;
        case 2: run(); break;
        case 3: timedOut(); break;
        case 4: readyReadSlot((*reinterpret_cast< const QHttpResponseHeader(*)>(_a[1]))); break;
        case 5: httpRequestFinished((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 6: httpDataReadProgress((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 7: httpDataSendProgress((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 8: httpRequestStarted((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: httpStateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: httpDone((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: readResponseHeader((*reinterpret_cast< const QHttpResponseHeader(*)>(_a[1]))); break;
        case 12: login(); break;
        case 13: execQuery(); break;
        case 14: getResults((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 15: followRedirect((*reinterpret_cast< const QHttpResponseHeader(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 16;
    }
    return _id;
}

// SIGNAL 0
void OpenMS::MascotRemoteQuery::done()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void OpenMS::MascotRemoteQuery::gotRedirect(const QHttpResponseHeader & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_END_MOC_NAMESPACE
