"""
fetch_all_solar_system.py
يحمل بيانات الموقع والسرعة والكتلة (GM) لأجسام المجموعة الشمسية الرئيسية
من NASA JPL Horizons ويكتبها في initial_all.csv.

الاستخدام: python3 fetch_all_solar_system.py 2026-01-01
المتطلبات: اتصال بالإنترنت فقط (بدون مكتبات خارجية)
"""

import sys
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timedelta

API = "https://ssd.jpl.nasa.gov/api/horizons.api"

# ============================================================
# قائمة الأجسام: (الاسم، معرّف Horizons، الأب، GM بوحدة km^3/s^2)
# GM للكواكب من DE440. GM للأقمار من جداول JPL.
# الكويكبات: نستخدم معرّفها الرقمي متبوعًا بـ ';'
# ============================================================
BODIES = [
    # ---- الشمس ----
    ("Sun",       "10",   "-",       132712440041.93938),

    # ---- الكواكب ----
    ("Mercury",   "199",  "Sun",     22031.868551),
    ("Venus",     "299",  "Sun",     324858.592000),
    ("Earth",     "399",  "Sun",     398600.435507),
    ("Mars",      "499",  "Sun",     42828.375816),
    ("Jupiter",   "599",  "Sun",     126712764.100000),
    ("Saturn",    "699",  "Sun",     37940584.841800),
    ("Uranus",    "799",  "Sun",     5794556.400000),
    ("Neptune",   "899",  "Sun",     6836527.100580),
    ("Pluto",     "999",  "Sun",     975.500000),

    # ---- القمر ----
    ("Moon",      "301",  "Earth",   4902.800118),

    # ---- أقمار المريخ ----
    ("Phobos",    "401",  "Mars",    0.0007112),
    ("Deimos",    "402",  "Mars",    0.0000985),

    # ---- أقمار المشتري (الجاليليون) ----
    ("Io",        "501",  "Jupiter", 5959.916),
    ("Europa",    "502",  "Jupiter", 3202.739),
    ("Ganymede",  "503",  "Jupiter", 9887.834),
    ("Callisto",  "504",  "Jupiter", 7179.289),

    # ---- أقمار زحل الرئيسية ----
    ("Mimas",     "601",  "Saturn",  2.5026),
    ("Enceladus", "602",  "Saturn",  7.2103),
    ("Tethys",    "603",  "Saturn",  41.2067),
    ("Dione",     "604",  "Saturn",  73.1146),
    ("Rhea",      "605",  "Saturn",  153.9417),
    ("Titan",     "606",  "Saturn",  8978.14),
    ("Iapetus",   "608",  "Saturn",  120.5038),

    # ---- أقمار أورانوس ----
    ("Ariel",     "701",  "Uranus",  83.5),
    ("Umbriel",   "702",  "Uranus",  85.1),
    ("Titania",   "703",  "Uranus",  228.2),
    ("Oberon",    "704",  "Uranus",  192.4),
    ("Miranda",   "705",  "Uranus",  4.4),

    # ---- أقمار نبتون ----
    ("Triton",    "801",  "Neptune", 1427.6),

    # ---- أقمار بلوتو ----
    ("Charon",    "901",  "Pluto",   105.88),

    # ============================================================
    # الكويكبات الكبيرة (معرّفها الرقمي + ';')
    # يمكنك إضافة المزيد هنا: ("اسم", "رقم;", "Sun", GM)
    # ============================================================
    ("Ceres",     "1;",   "Sun",     62.6284),
    ("Pallas",    "2;",   "Sun",     13.665),
    ("Vesta",     "4;",   "Sun",     17.288),
    ("Eros",      "433;", "Sun",     0.000446),
    ("Psyche",    "16;",  "Sun",     0.00153),
    ("Davida",    "511;", "Sun",     0.0043),
    ("Interamnia","704;", "Sun",     0.0021),
]


def ask_horizons(command, start, stop):
    """يطلب متجه الحالة لجسم واحد من Horizons ويعيد 6 أرقام."""
    params = {
        "format": "json",
        "COMMAND": "'" + command + "'",
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@0'",            # مركز ثقل المجموعة الشمسية
        "START_TIME": "'" + start + "'",
        "STOP_TIME": "'" + stop + "'",
        "STEP_SIZE": "'1d'",
        "VEC_TABLE": "'2'",
        "REF_SYSTEM": "'ICRF'",
        "REF_PLANE": "'ECLIPTIC'",
        "OUT_UNITS": "'AU-D'",
        "CSV_FORMAT": "'YES'",
    }
    url = API + "?" + urllib.parse.urlencode(params, quote_via=urllib.parse.quote)

    with urllib.request.urlopen(url, timeout=60) as response:
        text = json.loads(response.read().decode("utf-8"))["result"]

    if "$$SOE" not in text:
        raise RuntimeError(text[:400])

    first_line = text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()[0]
    fields = [f.strip() for f in first_line.split(",")]
    return [float(v) for v in fields[2:8]]


def gm_to_au_day2(gm_km3_s2):
    """تحويل GM من km^3/s^2 إلى AU^3/day^2."""
    AU_KM = 149597870.7
    DAY_S = 86400.0
    return gm_km3_s2 * (DAY_S ** 2) / (AU_KM ** 3)


def main():
    date_text = sys.argv[1] if len(sys.argv) > 1 else "2026-01-01"
    start_date = datetime.strptime(date_text, "%Y-%m-%d")
    start = start_date.strftime("%Y-%m-%d 00:00")
    stop = (start_date + timedelta(days=1)).strftime("%Y-%m-%d 00:00")

    rows = []
    for name, command, parent, gm in BODIES:
        print(f"Downloading {name} ...", end=" ", flush=True)
        try:
            state = ask_horizons(command, start, stop)
        except Exception as error:
            print("FAILED")
            print("   ", str(error)[:200])
            continue

        gm_au = gm_to_au_day2(gm)
        rows.append((name, parent, gm, gm_au, state))
        print("ok")
        time.sleep(0.3)

    with open("initial.csv", "w") as f:
        f.write("# epoch=" + start_date.strftime("%Y-%m-%dT00:00:00") + " (TDB)\n")
        f.write("name,parent,gm_km3_s2,gm_au3_d2,"
                "x_au,y_au,z_au,vx_aud,vy_aud,vz_aud\n")
        for name, parent, gm, gm_au, s in rows:
            f.write(f"{name},{parent},{gm!r},{gm_au:.15e},"
                    f"{s[0]:.15e},{s[1]:.15e},{s[2]:.15e},"
                    f"{s[3]:.15e},{s[4]:.15e},{s[5]:.15e}\n")

    print(f"\nSaved {len(rows)} bodies to initial_all.csv (epoch {date_text})")
    if len(rows) < len(BODIES):
        print("Some bodies failed. Run again, or remove them from BODIES.")


if __name__ == "__main__":
    main()