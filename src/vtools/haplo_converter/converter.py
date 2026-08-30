import os
import re
import json
import sqlite3
import logging
import sys

import requests
from difflib import get_close_matches

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def _app_dir() -> str:
    """Каталог, где лежит сам скрипт (или .exe при заморозке PyInstaller-ом).
    Используется как база для всех файлов данных, чтобы программа работала
    одинаково независимо от того, откуда её запустили (двойной клик,
    ярлык, другой рабочий каталог и т.п.)."""
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


APP_DIR = _app_dir()

YFULL_TREE_URL = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json"
LOCAL_TREE = os.path.join(APP_DIR, "current_tree.json")
DEFAULT_CSV = os.path.join(APP_DIR, "SNP-Index-Human.csv")
DEFAULT_DB_PATH = os.path.join(APP_DIR, "haplo_pro.db")

# Константы
ROOT_LABEL = "ROOT"
PATH_NOT_FOUND = "Путь не найден в дереве YFull"
MAX_RECURSION_DEPTH = 500

RE_YFULL = re.compile(r"^[A-Z]+-[A-Z0-9.]+$")
RE_ISOGG = re.compile(r"^[A-Z][A-Z0-9]+$")
RE_ALNUM = re.compile(r"^[0-9A-Z.]+$")

ISOGG_TO_YFULL = {

    # ──────────────── HAPLOGROUP A (корень) ────────────────
    "A":            "A-Y119988",    # самый корень Y-дерева
    "A0":           "A-L1085",
    "A00":          "A-L1149",
    "A0a":          "A-P108",
    "A0b":          "A-V148",
    "A1":           "A-L985",
    "A1a":          "A-M31",
    "A1b":          "A-P108",
    "A1b1":         "A-M32",
    "A2":           "A-V50",
    "A3":           "A-M144",

    # ──────────────── HAPLOGROUP B ────────────────
    "B":            "B-M60",
    "B2":           "B-M182",
    "B2a":          "B-M150",
    "B2b":          "B-M112",

    # ──────────────── HAPLOGROUP C ────────────────
    "C":            "C-M130",
    "C1":           "C-CTS11043",
    "C1a":          "C-CTS2657",
    "C1a1":         "C-M8",
    "C1b":          "C-F3393",
    "C1b1":         "C-K281",
    "C1b2":         "C-M38",
    "C1b2a":        "C-M208",
    "C1b2b":        "C-M347",
    "C2":           "C-M217",
    "C2a":          "C-L1373",
    "C2b":          "C-M48",
    "C2b1":         "C-M86",
    "C2b1a":        "C-F1396",
    "C2c":          "C-M407",
    "C3":           "C-M77",

    # ──────────────── HAPLOGROUP D ────────────────
    "D":            "D-M174",
    "D1":           "D-CTS3946",
    "D1a":          "D-M15",
    "D1a1":         "D-P99",
    "D1a2":         "D-M55",
    "D1a2a":        "D-M116.1",
    "D1b":          "D-M64.1",

    # ──────────────── HAPLOGROUP E ────────────────
    "E":            "E-M5479",      # ← ваш исходный
    "E1":           "E-P147",
    "E1a":          "E-M33",
    "E1a1":         "E-M44",
    "E1b":          "E-P177",
    "E1b1":         "E-P2",
    "E1b1a":        "E-V38",
    "E1b1a1":       "E-M2",
    "E1b1a1a":      "E-M58",
    "E1b1a1b":      "E-M116.2",
    "E1b1a1c":      "E-M149",
    "E1b1a1d":      "E-M154",
    "E1b1a1e":      "E-M191",
    "E1b1a1e1":     "E-U186",
    "E1b1a1f":      "E-Z1123",
    "E1b1a1g":      "E-U175",
    "E1b1a1g1":     "E-U209",
    "E1b1a1g2":     "E-U290",
    "E1b1b":        "E-M215",
    "E1b1b1":       "E-M35",
    "E1b1b1a":      "E-M78",
    "E1b1b1a1":     "E-V12",
    "E1b1b1a1a":    "E-V32",
    "E1b1b1a1b":    "E-V22",
    "E1b1b1a2":     "E-V13",
    "E1b1b1a3":     "E-V65",
    "E1b1b1b":      "E-M81",
    "E1b1b1b1":     "E-M107",
    "E1b1b1b2":     "E-M183",
    "E1b1b1c":      "E-M123",
    "E1b1b1c1":     "E-M34",
    "E1b1b1d":      "E-M281",
    "E1b1b1e":      "E-V6",
    "E2":           "E-M75",
    "E2a":          "E-M41",
    "E2b":          "E-M90",

    # ──────────────── HAPLOGROUP F (ствол) ────────────────
    "F":            "F-M89",
    "F1":           "F-P91",
    "F2":           "F-M427",
    "F3":           "F-M481",

    # ──────────────── HAPLOGROUP G ────────────────
    "G":            "G-M201",       # ← ваш исходный
    "G1":           "G-M285",
    "G1a":          "G-P20",
    "G1b":          "G-L1323",
    "G2":           "G-P287",
    "G2a":          "G-P15",
    "G2a1":         "G-FGC7535",
    "G2a2":         "G-L30",
    "G2a2a":        "G-PF3147",
    "G2a2b":        "G-L43",
    "G2a2b1":       "G-L140",
    "G2a2b2":       "G-PF3267",
    "G2a2b2a":      "G-L497",
    "G2a2b2a1":     "G-L42",
    "G2a2b2a1a":    "G-U1",
    "G2a2b2a1b":    "G-CTS9737",
    "G2b":          "G-M286",
    "G2b1":         "G-M377",
    "G2b2":         "G-M283",

    # ──────────────── HAPLOGROUP H ────────────────
    "H":            "H-L901",
    "H1":           "H-M69",
    "H1a":          "H-M52",
    "H1a1":         "H-M82",
    "H2":           "H-P96",
    "H3":           "H-Z5857",

    # ──────────────── HAPLOGROUP I ────────────────
    "I":            "I-M170",       # ← ваш исходный
    "I1":           "I-DF29",       # ← ваш исходный
    "I1a":          "I-M253",
    "I1a1":         "I-DF29",       # (DF29 = определяющий SNP I1 в YFull)
    "I1a1a":        "I-Z60",
    "I1a1a1":       "I-Z63",
    "I1a1a2":       "I-Z58",
    "I1a1b":        "I-Z140",
    "I1a2":         "I-S244",
    "I1a3":         "I-Z74",
    "I2":           "I-S31",
    "I2a":          "I-L460",
    "I2a1":         "I-P37.2",
    "I2a1a":        "I-CTS10228",
    "I2a1a1":       "I-L147.2",
    "I2a1a2":       "I-CTS5966",
    "I2a1b":        "I-M423",
    "I2a1b1":       "I-CTS10936",
    "I2a1b1a":      "I-Y3120",
    "I2a1b2":       "I-CTS1977",
    "I2a2":         "I-M436",
    "I2a2a":        "I-M223",
    "I2a2a1":       "I-L701",
    "I2a2a1a":      "I-CTS6433",
    "I2a2a1b":      "I-CTS616",
    "I2a2b":        "I-L38",
    "I2b":          "I-L416",
    "I2c":          "I-L596",

    # ──────────────── HAPLOGROUP J ────────────────
    "J":            "J-M304",       # ← ваш исходный
    "J1":           "J-M267",
    "J1a":          "J-P209",
    "J1a1":         "J-L255",
    "J1a1a":        "J-L147",
    "J1a2":         "J-Z1828",
    "J1a2a":        "J-L620",
    "J1a2a1":       "J-Z2223",
    "J1a2a1a":      "J-Y6304",
    "J1a2b":        "J-Y5575",
    "J2":           "J-M172",
    "J2a":          "J-M410",
    "J2a1":         "J-L26",
    "J2a1a":        "J-DYS445=6",
    "J2a1b":        "J-M67",
    "J2a1b1":       "J-M92",
    "J2a1h":        "J-L24",
    "J2a1h1":       "J-L25",
    "J2b":          "J-M12",
    "J2b1":         "J-M102",
    "J2b1a":        "J-Z1296",
    "J2b2":         "J-M241",

    # ──────────────── HAPLOGROUP K (ствол) ────────────────
    "K":            "K-M9",
    "K2":           "K-M526",
    "K2a":          "K-M2335",
    "K2b":          "K-P331",
    "K2b1":         "K-P397",
    "K2c":          "K-P261",
    "K2d":          "K-P402",
    "K2e":          "K-M147",

    # ──────────────── HAPLOGROUP L ────────────────
    "L":            "L-M20",
    "L1":           "L-M76",
    "L1a":          "L-M27",
    "L1b":          "L-M317",
    "L1c":          "L-M357",
    "L2":           "L-L595",

    # ──────────────── HAPLOGROUP M ────────────────
    "M":            "M-P256",
    "M1":           "M-M4",
    "M2":           "M-M353",

    # ──────────────── HAPLOGROUP N ────────────────
    "N":            "N-M231",       # ← ваш исходный
    "N1":           "N-CTS2929",
    "N1a":          "N-L708",
    "N1a1":         "N-M128",
    "N1a2":         "N-P43",
    "N1a2a":        "N-P63",
    "N1b":          "N-P63",
    "N1c":          "N-M46",
    "N1c1":         "N-M178",
    "N1c1a":        "N-L1026",
    "N1c1a1":       "N-Z4908",
    "N1c1a1a":      "N-CTS3103",
    "N1c1a1a1":     "N-Y9022",
    "N1c1a1a2":     "N-VL29",

    # ──────────────── HAPLOGROUP O ────────────────
    "O":            "O-M175",
    "O1":           "O-CTS2780",
    "O1a":          "O-M119",
    "O1a1":         "O-P203",
    "O1b":          "O-P31",
    "O1b1":         "O-M95",
    "O1b1a":        "O-M88",
    "O1b2":         "O-P49",
    "O2":           "O-M122",
    "O2a":          "O-M324",
    "O2a1":         "O-M159",
    "O2a2":         "O-M188",
    "O2a2a":        "O-M7",
    "O2a2b":        "O-M134",
    "O2a2b1":       "O-M117",
    "O2a2b1a":      "O-F444",
    "O2b":          "O-P49",

    # ──────────────── HAPLOGROUP P ────────────────
    "P":            "P-P295",
    "P1":           "P-M45",
    "P1a":          "P-M74",
    "P2":           "P-P226",

    # ──────────────── HAPLOGROUP Q ────────────────
    "Q":            "Q-M242",
    "Q1":           "Q-F903",
    "Q1a":          "Q-F746",
    "Q1a1":         "Q-M120",
    "Q1a2":         "Q-M25",
    "Q1a2a":        "Q-L712",
    "Q1a2a1":       "Q-M346",
    "Q1a2a1a":      "Q-L54",
    "Q1a2a1a1":     "Q-M3",
    "Q1a2a1a2":     "Q-Z780",
    "Q1a3":         "Q-L330",
    "Q1b":          "Q-M378",
    "Q1b1":         "Q-L245",
    "Q1b1a":        "Q-M323",

    # ──────────────── HAPLOGROUP R ────────────────
    "R":            "R-M207",
    "R1":           "R-M173",
    "R1a":          "R-Y482",        # ← ваш исходный
    "R1a1":         "R-M459",        # ← ваш исходный
    "R1a1a":        "R-M17",         # главная ветвь R1a (M17/M198)
    "R1a1a1":       "R-M417",        # корень всех субкладов R1a
    "R1a1a1a":      "R-L664",
    "R1a1a1b":      "R-Z645",
    "R1a1a1b1":     "R-Z283",
    "R1a1a1b1a":    "R-Z282",
    "R1a1a1b1a1":   "R-Z284",
    "R1a1a1b1a2":   "R-CTS1211",
    "R1a1a1b1a2a":  "R-Z91",
    "R1a1a1b1a2b":  "R-CTS3402",
    "R1a1a1b1a3":   "R-Z280",
    "R1a1a1b1a3a":  "R-CTS1954",
    "R1a1a1b1b":    "R-Z93",
    "R1a1a1b1b1":   "R-Z94",
    "R1a1a1b1b1a":  "R-Y57",
    "R1a1a1b1b2":   "R-Y934",
    "R1a1a1b1b3":   "R-YP4800",
    "R1a1a1b2":     "R-Z2",
    "R1a1a1b2a":    "R-Z2124",
    "R1a1a2":       "R-Y19144",
    "R1b":          "R-L754",        # ← ваш исходный
    "R1b1":         "R-L761",        # ← ваш исходный
    "R1b1a":        "R-L389",        # ← ваш исходный
    "R1b1a1":       "R-P297",
    "R1b1a1a":      "R-M73",
    "R1b1a1b":      "R-M269",
    "R1b1a1b1":     "R-L23",
    "R1b1a1b1a":    "R-L51",
    "R1b1a1b1a1":   "R-P310",
    "R1b1a1b1a1a":  "R-U106",
    "R1b1a1b1a1a1": "R-Z381",
    "R1b1a1b1a1a2": "R-Z156",
    "R1b1a1b1a1b":  "R-DF27",
    "R1b1a1b1a1b1": "R-Z195",
    "R1b1a1b1a2":   "R-L21",
    "R1b1a1b1a2a":  "R-DF13",
    "R1b1a1b1a2a1": "R-Z251",
    "R1b1a1b1a2b":  "R-DF21",
    "R1b1a1b1b":    "R-Z2103",
    "R1b1a1b1b1":   "R-Z2105",
    "R1b1a1b1b2":   "R-CTS7556",
    "R1b1a1b1b3":   "R-PF7562",
    "R1b1a2":       "R-V88",
    "R1b1a2a":      "R-M18",
    "R1b1b":        "R-PH155",
    "R1b2":         "R-PH200",
    "R2":           "R-M479",
    "R2a":          "R-M124",

    # ──────────────── HAPLOGROUP S ────────────────
    "S":            "S-M230",
    "S1":           "S-P405",
    "S1a":          "S-M254",

    # ──────────────── HAPLOGROUP T ────────────────
    "T":            "T-M184",       # ← ваш исходный
    "T1":           "T-M272",
    "T1a":          "T-M70",
    "T1a1":         "T-L162",
    "T1a1a":        "T-L208",
    "T1a1a1":       "T-CTS2214",
    "T1a2":         "T-P321",

}

YFULL_TO_ISOGG = {v.upper(): k.upper() for k, v in ISOGG_TO_YFULL.items()}

# Case-insensitive lookup index: dict keys above use mixed case (e.g. "R1a"),
# but user input and CSV/JSON data arrive in varying case. Build an
# upper-cased lookup table once so exact matches are not silently missed.
ISOGG_TO_YFULL_UPPER = {k.upper(): v for k, v in ISOGG_TO_YFULL.items()}

# ================= DATABASE LAYER =================

class HaploDB:
    def __init__(self, db_path=None):
        if db_path is None:
            db_path = DEFAULT_DB_PATH
        db_dir = os.path.dirname(db_path)
        if db_dir and not os.path.isdir(db_dir):
            os.makedirs(db_dir, exist_ok=True)
        try:
            self.conn = sqlite3.connect(db_path, check_same_thread=False)
        except sqlite3.OperationalError as e:
            raise RuntimeError(
                f"Не удалось открыть базу данных по пути:\n{db_path}\n"
                f"Причина: {e}\n"
                "Проверьте, что папка существует и у программы есть права на запись в неё "
                "(например, запустите не из системной/защищённой папки типа Program Files)."
            ) from e
        self._setup()

    def _setup(self):
        cursor = self.conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS snp_data (
                name   TEXT PRIMARY KEY,
                hg     TEXT,
                path   TEXT,
                source TEXT,
                extra  TEXT
            )
        """)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_snp_name ON snp_data(name)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_snp_hg ON snp_data(hg)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_snp_path ON snp_data(path)")
        self.conn.commit()

    def clear(self):
        self.conn.execute("DELETE FROM snp_data")
        self.conn.commit()

    def get_stats(self):
        return self.conn.execute("SELECT COUNT(*) FROM snp_data").fetchone()[0]


# ================= LOGIC LAYER =================

class Converter:
    def __init__(self, db: HaploDB):
        self.db = db

    def normalize(self, text):
        if not text:
            return ""
        return re.sub(r"[\s\-]", "", str(text)).upper()

    def download_yfull_tree(self, url=YFULL_TREE_URL, local_path=LOCAL_TREE):
        try:
            logger.info(f"Downloading YFull tree from {url}")
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            with open(local_path, "w", encoding="utf-8") as f:
                f.write(resp.text)
            logger.info(f"YFull tree saved here: {local_path}")
            return local_path
        except Exception as e:
            logger.error(f"Downloading YFull tree failed: {e}")
            if os.path.exists(local_path):
                logger.info(f"Use a local cash: {local_path}")
                return local_path
            raise RuntimeError("Unable to retrieve YFull tree either over the network or from the local cache")

    def import_yfull_online(self):
        path = self.download_yfull_tree()
        self.import_yfull_json(path)

    def find_best_path_by_prefix(self, yfull_hg):
        """
        Ищет лучший полный путь в YFull-данных:
        - сначала точный hg,
        - потом вхождение в path,
        - потом поднимается к более короткому предку.
        """
        target = str(yfull_hg).strip().upper()
        if not target:
            return PATH_NOT_FOUND

        candidates = [target]

        # если R-L51 -> потом R
        if "-" in target:
            head = target.split("-", 1)[0]
            if head and head not in candidates:
                candidates.append(head)

        for cand in candidates:
            row = self.db.conn.execute(
                """
                SELECT path FROM snp_data
                WHERE source = 'YFull_JSON'
                  AND path IS NOT NULL
                  AND path != ''
                  AND (
                        UPPER(hg) = ?
                     OR UPPER(path) LIKE ?
                     OR UPPER(path) LIKE ?
                     OR UPPER(path) LIKE ?
                  )
                ORDER BY LENGTH(path) DESC
                LIMIT 1
                """,
                (
                    cand,
                    f"%/{cand}/%",
                    f"%/{cand}",
                    f"{ROOT_LABEL}/{cand}%"
                )
            ).fetchone()

            if row and row[0]:
                return row[0]

        return PATH_NOT_FOUND
    
    def get_yfull_equivalent(self, isogg_hg):
        if not isogg_hg:
            return ""
        isogg = str(isogg_hg).strip().upper()
        if isogg in ISOGG_TO_YFULL_UPPER:
            return ISOGG_TO_YFULL_UPPER[isogg]

        parent = isogg
        while len(parent) > 1:
            parent = parent[:-1]
            if parent in ISOGG_TO_YFULL_UPPER:
                return ISOGG_TO_YFULL_UPPER[parent]
        return ""

    def format_slash_path(self, path_value):
        if not path_value or path_value == PATH_NOT_FOUND:
            return PATH_NOT_FOUND

        value = str(path_value).strip()

        if "/" not in value:
            found = self.find_best_path_by_prefix(value)
            if found != PATH_NOT_FOUND:
                value = found
            else:
                return value

        parts = [p.strip() for p in value.split("/") if p.strip()]
        parts = [p for p in parts if p.upper() != ROOT_LABEL]

        if not parts:
            return PATH_NOT_FOUND

        return "/" + "/".join(parts)

    def search_by_haplogroup_fast(self, query_hg):
        q = str(query_hg).strip().upper()
        if not q:
            return None

        # 1. Точная запись в базе
        row = self.db.conn.execute(
            """
            SELECT hg, path, source FROM snp_data
            WHERE UPPER(hg) = ?
            LIMIT 1
            """,
            (q,)
        ).fetchone()

        if row:
            hg, path, source = row
            real_path = path
            if not real_path or real_path in ("CSV_ISOGG_INDEX", PATH_NOT_FOUND):
                real_path = self.resolve_csv_hg_path(hg)

            return {
                "query": q,
                "hg_yfull": self.get_yfull_equivalent(q) or hg,
                "path": self.format_slash_path(real_path),
                "source": source
            }

        # 2. Через словарь ISOGG -> YFull
        yfull = self.get_yfull_equivalent(q)
        if yfull:
            real_path = self.find_best_path_by_prefix(yfull)
            return {
                "query": q,
                "hg_yfull": yfull,
                "path": self.format_slash_path(real_path),
                "source": "ISOGG_TO_YFULL"
            }

        return None
    
    def resolve_csv_hg_path(self, isogg_hg):
        if not isogg_hg:
            return PATH_NOT_FOUND

        isogg = str(isogg_hg).strip().upper()

        # 1. Точное совпадение
        yfull_hg = ISOGG_TO_YFULL_UPPER.get(isogg)
        if yfull_hg:
            found = self.find_best_path_by_prefix(yfull_hg)
            if found != PATH_NOT_FOUND:
                return found

        # 2. Ближайший родитель справа налево
        parent = isogg
        while len(parent) > 1:
            parent = parent[:-1]
            yfull_hg = ISOGG_TO_YFULL_UPPER.get(parent)
            if yfull_hg:
                found = self.find_best_path_by_prefix(yfull_hg)
                if found != PATH_NOT_FOUND:
                    return found

        return PATH_NOT_FOUND

    def format_full_path(self, path_value):
        if not path_value or path_value == PATH_NOT_FOUND:
            return PATH_NOT_FOUND

        value = str(path_value).strip()

        if "/" not in value:
            found = self.find_best_path_by_prefix(value)
            if found != PATH_NOT_FOUND:
                value = found
            else:
                return value

        parts = [p.strip() for p in value.split("/") if p.strip()]
        parts = [p for p in parts if p.upper() != ROOT_LABEL]

        if not parts:
            return PATH_NOT_FOUND

        return " -> ".join(parts)

    def get_adam_path(self, hg_or_path):
        if not hg_or_path:
            return PATH_NOT_FOUND

        value = str(hg_or_path).strip()
        if not value:
            return PATH_NOT_FOUND

        if "/" in value:
            parts = [p for p in value.split("/") if p]
        else:
            parts = [value]

        # убираем ROOT из начала
        if parts and parts[0].upper() == ROOT_LABEL:
            parts = parts[1:]

        if not parts:
            return PATH_NOT_FOUND

        return " -> ".join(parts)

    def import_yfull_json(self, path):
        with open(path, encoding="utf-8") as f:
            tree = json.load(f)

        batch = []
        stack = [(tree, [])]
        while stack:
            node, path_list = stack.pop()
            hg = str(node.get("id", ROOT_LABEL)).strip().upper()
            new_path = path_list + [hg]
            path_str = "/".join(new_path)

            snps = node.get("snps", "")
            for s in re.split(r"[,/]", snps):
                s_norm = self.normalize(s)
                if s_norm:
                    batch.append((s_norm, hg, path_str, "YFull_JSON", "{}"))

            for c in node.get("children", []):
                stack.append((c, new_path))

        self.db.conn.executemany("INSERT OR REPLACE INTO snp_data VALUES (?,?,?,?,?)", batch)
        self.db.conn.commit()

    def download_default_csv(self, csv_url=None, local_path=DEFAULT_CSV):
        if csv_url is None:
            if os.path.exists(local_path):
                return local_path
            raise RuntimeError("Default CSV file not found")

        try:
            logger.info(f"Downloading CSV from {csv_url}")
            resp = requests.get(csv_url, timeout=30)
            resp.raise_for_status()
            with open(local_path, "wb") as f:
                f.write(resp.content)
            logger.info(f"CSV saved here: {local_path}")
            return local_path
        except Exception as e:
            logger.error(f"Downloading CSV failed: {e}")
            if os.path.exists(local_path):
                logger.info(f"Use a local CSV: {local_path}")
                return local_path
            raise RuntimeError("Unable to retrieve CSV either over the network or from the local cache")

    def looks_like_haplogroup(self, text):
        """
        True только для запросов вида:
        R1B1A1B1A1A2C1A1D5A, J2, N1C1, I2A1 и т.п.
        Но НЕ для SNP вроде A100, M269, Z280, CTS1211, FGC8628.
        """
        if not text:
            return False

        q = str(text).strip().upper()

        if re.fullmatch(r"[A-Z]+-[A-Z0-9.]+", q):
            return True

        if re.fullmatch(r"[A-Z][0-9A-Z]+", q):
            snp_like_prefixes = (
                "A", "M", "Z", "CTS", "FGC", "L", "P", "V", "Y", "BY", "FT",
                "PF", "SK", "PH", "S", "U", "F", "YP", "ZZ", "MF", "DF"
            )
            for pref in snp_like_prefixes:
                if q.startswith(pref):
                    if q == pref:
                        continue
                    tail = q[len(pref):]
                    if tail and re.fullmatch(r"[0-9A-Z.]+", tail):
                        return False

            return True

        return False  

    def search_by_haplogroup_fast(self, query_hg):
        q = str(query_hg).strip().upper()
        if not q:
            return None

        if not self.looks_like_haplogroup(q):
            return None

        row = self.db.conn.execute(
            """
            SELECT hg, path, source FROM snp_data
            WHERE UPPER(hg) = ?
            LIMIT 1
            """,
            (q,)
        ).fetchone()

        if row:
            hg, path, source = row
            real_path = path
            if not real_path or real_path in ("CSV_ISOGG_INDEX", PATH_NOT_FOUND):
                real_path = self.resolve_csv_hg_path(hg)

            return {
                "query": q,
                "display_hg": q,
                "hg_yfull": self.get_yfull_equivalent(q) or hg,
                "path": self.format_slash_path(real_path),
                "source": source
            }

        yfull = self.get_yfull_equivalent(q)
        if yfull:
            real_path = self.find_best_path_by_prefix(yfull)
            return {
                "query": q,
                "display_hg": q,
                "hg_yfull": yfull,
                "path": self.format_slash_path(real_path),
                "source": "ISOGG_TO_YFULL"
            }

        return None      

    def import_csv(self, path):
        import csv

        def clean_col(s):
            return str(s).replace("\ufeff", "").strip().strip('"').strip()

        with open(path, "r", encoding="utf-8-sig", newline="") as f:
            rows = list(csv.reader(f))

        if not rows:
            raise ValueError("CSV is empty.")

        header_index = None
        for i, row in enumerate(rows[:100]):
            cleaned = [clean_col(x).upper() for x in row]
            if "SUBGROUP NAME" in cleaned and "NAME" in cleaned:
                header_index = i
                break

        if header_index is None:
            preview = "\n".join([",".join(r[:5]) for r in rows[:10]])
            raise ValueError(
                "Unable to find a header row with columns 'Subgroup Name' and 'Name'.\n\n"
                f"First lines of file:\n{preview}"
            )

        header = [clean_col(x) for x in rows[header_index]]
        header_upper = [x.upper() for x in header]

        try:
            hg_idx = header_upper.index("SUBGROUP NAME")
            snp_idx = header_upper.index("NAME")
        except ValueError:
            raise ValueError(f"Header row found, but columns not recognized: {header}")

        batch = []
        for row in rows[header_index + 1:]:
            if not row or not any(str(x).strip() for x in row):
                continue

            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            elif len(row) > len(header):
                row = row[:len(header)]

            snp_raw = row[snp_idx]
            hg_raw = row[hg_idx]

            name = self.normalize(snp_raw)
            hg = str(hg_raw).strip().upper()

            if not name or not hg:
                continue

            extra_dict = {}
            for idx, value in enumerate(row):
                if idx in (snp_idx, hg_idx):
                    continue
                extra_dict[header[idx]] = value

            extra = json.dumps(extra_dict, ensure_ascii=False)
            resolved_path = self.resolve_csv_hg_path(hg)
            batch.append((name, hg, resolved_path, os.path.basename(path), extra))

        if not batch:
            raise ValueError("CSV was found, but no valid records could be extracted.")

        self.db.conn.executemany(
            "INSERT OR REPLACE INTO snp_data VALUES (?,?,?,?,?)",
            batch
        )
        self.db.conn.commit()

    def detect_hg_notation(self, value):
        """
        Возвращает:
        - 'YFULL' для R-L51, I-M170, N-M231
        - 'ISOGG' для R1B1A1B1A1A2...
        - None если не распознано
        """
        if not value:
            return None

        v = str(value).strip().upper()

        if RE_YFULL.fullmatch(v):
            return "YFULL"

        if RE_ISOGG.fullmatch(v):
            # простые SNP не считаем ISOGG
            snp_like_prefixes = (
                "A", "M", "Z", "CTS", "FGC", "L", "P", "V", "Y", "BY", "FT",
                "PF", "SK", "PH", "S", "U", "F", "YP", "ZZ", "MF", "DF"
            )
            for pref in snp_like_prefixes:
                if v.startswith(pref):
                    if v == pref:
                        continue
                    tail = v[len(pref):]
                    if tail and RE_ALNUM.fullmatch(tail):
                        return None
            return "ISOGG"

        return None

    def get_yfull_equivalent(self, isogg_hg):
        if not isogg_hg:
            return ""
        isogg = str(isogg_hg).strip().upper()

        if isogg in ISOGG_TO_YFULL_UPPER:
            return ISOGG_TO_YFULL_UPPER[isogg]

        parent = isogg
        while len(parent) > 1:
            parent = parent[:-1]
            if parent in ISOGG_TO_YFULL_UPPER:
                return ISOGG_TO_YFULL_UPPER[parent]

        return ""

    def get_isogg_equivalent(self, yfull_hg):
        if not yfull_hg:
            return ""
        yfull = str(yfull_hg).strip().upper()

        if yfull in YFULL_TO_ISOGG:
            return YFULL_TO_ISOGG[yfull]

        # если точного совпадения нет, пробуем подняться к букве-основе
        if "-" in yfull:
            head = yfull.split("-", 1)[0]
            if head in YFULL_TO_ISOGG:
                return YFULL_TO_ISOGG[head]

        return ""

    def convert_haplogroup_value(self, value):
        """
        Возвращает:
        {
            'input': исходное,
            'detected': ISOGG/YFULL/None,
            'converted': результат или '',
            'direction': 'ISOGG->YFULL' / 'YFULL->ISOGG' / ''
        }
        """
        raw = str(value).strip()
        up = raw.upper()

        detected = self.detect_hg_notation(up)
        if detected == "ISOGG":
            converted = self.get_yfull_equivalent(up)
            return {
                "input": raw,
                "detected": "ISOGG",
                "converted": converted,
                "direction": "ISOGG->YFULL" if converted else ""
            }

        if detected == "YFULL":
            converted = self.get_isogg_equivalent(up)
            return {
                "input": raw,
                "detected": "YFULL",
                "converted": converted,
                "direction": "YFULL->ISOGG" if converted else ""
            }

        return {
            "input": raw,
            "detected": None,
            "converted": "",
            "direction": ""
        }

    def convert_file(self, input_path, progress_callback=None):
        import csv

        output_path = self.make_converted_filename(input_path)
        results = []
        cache = {}

        # Сначала считаем примерный объём
        with open(input_path, "r", encoding="utf-8-sig", newline="") as f:
            sample = f.read(4096)
            f.seek(0)

            try:
                dialect = csv.Sniffer().sniff(sample)
            except Exception:
                dialect = csv.excel

            reader = csv.reader(f, dialect)
            total_cells = 0
            for row in reader:
                if isinstance(row, list):
                    total_cells += len(row)

        if total_cells == 0:
            raise ValueError("Файл пустой или не содержит данных.")

        done = 0

        with open(input_path, "r", encoding="utf-8-sig", newline="") as f:
            reader = csv.reader(f, dialect)

            for row_idx, row in enumerate(reader, start=1):
                if not isinstance(row, list) or not row:
                    continue

                for col_idx, cell in enumerate(row, start=1):
                    cell_key = str(cell).strip().upper()

                    if cell_key in cache:
                        conv = cache[cell_key]
                    else:
                        conv = self.convert_haplogroup_value(cell)
                        cache[cell_key] = conv

                    if conv["detected"] and conv["converted"]:
                        results.append({
                            "row": row_idx,
                            "column": col_idx,
                            "original": conv["input"],
                            "detected_format": conv["detected"],
                            "direction": conv["direction"],
                            "converted": conv["converted"]
                        })

                    done += 1
                    if progress_callback and (done % 500 == 0 or done == total_cells):
                        progress_callback(done, total_cells, f"Обработка... {done}/{total_cells}")

        if not results:
            raise ValueError("В файле не найдено ни одной распознанной гаплогруппы для конвертации.")

        with open(output_path, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "row",
                    "column",
                    "original",
                    "detected_format",
                    "direction",
                    "converted"
                ]
            )
            writer.writeheader()
            writer.writerows(results)

        return output_path, len(results)
    
    def make_converted_filename(self, input_path):
        base, ext = os.path.splitext(input_path)
        if not ext:
            ext = ".csv"
        return f"{base}_converted{ext if ext.lower() in ('.csv', '.txt') else '.csv'}"


# ================= CLI LAYER =================

class Cli:
    """Консольное меню, полностью повторяющее возможности прежнего GUI:
    загрузка YFull/JSON/CSV, конвертация файла, очистка БД, поиск."""

    def __init__(self):
        self.db = HaploDB()
        self.conv = Converter(self.db)

        # Автоматическая инициализация из сети, если база пуста (как в GUI)
        try:
            if self.db.get_stats() == 0:
                print("База данных пуста — пытаюсь скачать дерево YFull из сети...")
                self.conv.import_yfull_online()
                print(f"✅ Дерево YFull загружено. Записей в базе: {self.db.get_stats()}")
        except Exception as e:
            logger.error(f"Автоинициализация не удалась: {e}")
            print(f"⚠️ Не удалось автоматически загрузить дерево YFull: {e}")
            print("   Вы можете сделать это вручную через пункт меню.")

    # ---------- вспомогательные ----------
    def update_status(self):
        count = self.db.get_stats()
        print(f"📊 Записей в базе: {count}")

    def prompt_path(self, prompt_text: str) -> str | None:
        raw = input(prompt_text).strip().strip('"')
        if not raw:
            return None
        path = os.path.abspath(raw)
        if not os.path.isfile(path):
            print(f"❌ Файл не найден: {path}")
            return None
        return path

    # ---------- пункты меню ----------
    def load_yfull_online(self):
        try:
            print("🌐 Скачиваю дерево YFull из сети...")
            self.conv.import_yfull_online()
            self.update_status()
            print("✅ Дерево YFull скачано и добавлено в базу!")
        except Exception as e:
            print(f"❌ Ошибка: {e}")

    def load_file(self, mode):
        label = "JSON" if mode == "json" else "CSV"
        path = self.prompt_path(f"Введите путь к {label}-файлу: ")
        if not path:
            return
        try:
            if mode == "json":
                self.conv.import_yfull_json(path)
            else:
                self.conv.import_csv(path)
            self.update_status()
            print("✅ Данные добавлены в базу!")
        except Exception as e:
            print(f"❌ Ошибка: {e}")

    def clear_db(self):
        answer = input("Удалить все данные из базы? [y/N]: ").strip().lower()
        if answer in ("y", "yes"):
            self.db.clear()
            self.update_status()
            print("🗑 База очищена.")
        else:
            print("Отменено.")

    def convert_file_cli(self):
        path = self.prompt_path("Введите путь к файлу для конвертации (.txt/.csv): ")
        if not path:
            return

        last_percent = {"value": -1}

        def progress_callback(done, total, status_text):
            percent = (done / total) * 100 if total else 0
            # печатаем не чаще, чем на каждый процент, чтобы не спамить терминал
            if int(percent) != last_percent["value"]:
                last_percent["value"] = int(percent)
                bar_len = 30
                filled = int(bar_len * percent / 100)
                bar = "█" * filled + "░" * (bar_len - filled)
                sys.stdout.write(f"\r[{bar}] {percent:5.1f}%  {status_text}")
                sys.stdout.flush()

        try:
            print("🔄 Конвертация начата...")
            output_path, count = self.conv.convert_file(path, progress_callback=progress_callback)
            print()  # перевод строки после прогресс-бара
            print(f"✅ Готово. Конвертировано записей: {count}")
            print(f"📁 Файл сохранён: {output_path}")
        except Exception as e:
            print()
            print(f"❌ Ошибка конвертации: {e}")

    def search(self):
        raw_query = input("Введите SNP, гаплогруппу или подстроку для поиска: ").strip()
        if not raw_query:
            return
        try:
            query = self.conv.normalize(raw_query)
            if not query:
                return

            print()
            cursor = self.db.conn.cursor()

            # 1. Быстрый режим ТОЛЬКО для гаплогрупп
            fast_result = self.conv.search_by_haplogroup_fast(raw_query)
            if fast_result:
                print(f"🧬 Search: {fast_result['query']}")
                print(f"   Haplogroup: {fast_result['display_hg']}")
                print(f"   Path: {fast_result['path']}")
                print(f"   Source: {fast_result['source']}")
                if fast_result.get("hg_yfull"):
                    print(f"   YFull-eqv.: {fast_result['hg_yfull']}")
                print("\n" + "-" * 50 + "\n")
                return

            # 2. Обычный поиск SNP / мутаций
            cursor.execute("SELECT * FROM snp_data WHERE UPPER(name) = ?", (query,))
            rows = cursor.fetchall()

            # 3. Поиск по hg, path, source
            if not rows:
                like_q = f"%{raw_query.upper()}%"
                cursor.execute("""
                    SELECT * FROM snp_data
                    WHERE UPPER(name) LIKE ?
                       OR UPPER(hg) LIKE ?
                       OR UPPER(path) LIKE ?
                       OR UPPER(source) LIKE ?
                    LIMIT 100
                """, (like_q, like_q, like_q, like_q))
                rows = cursor.fetchall()

            # 4. Нечёткий поиск только по name
            if not rows:
                cursor.execute("SELECT name FROM snp_data LIMIT 50000")
                all_names = [r[0] for r in cursor.fetchall()]
                close = get_close_matches(query, all_names, n=10, cutoff=0.75)
                if close:
                    placeholder = ",".join("?" for _ in close)
                    cursor.execute(f"SELECT * FROM snp_data WHERE name IN ({placeholder})", close)
                    rows = cursor.fetchall()

            if not rows:
                print(f"❌ Ничего не найдено по запросу '{raw_query}'.")
                return

            for row in rows:
                if len(row) == 5:
                    name, hg, path, src, extra = row
                else:
                    _, name, hg, path, src, extra = row

                real_path = path
                if not real_path or real_path in ("CSV_ISOGG_INDEX", PATH_NOT_FOUND):
                    real_path = self.conv.resolve_csv_hg_path(hg)

                full_path = self.conv.format_slash_path(real_path)

                print(f"🧬 Search: {name}")
                print(f"   Гаплогруппа: {hg}")
                print(f"   Путь: {full_path}")
                print(f"   Источник: {src}")

                yfull_equiv = self.conv.get_yfull_equivalent(hg)
                if yfull_equiv:
                    print(f"   YFull-эквивалент: {yfull_equiv}")

                try:
                    ext_data = json.loads(extra) if extra else {}
                    if isinstance(ext_data, dict) and ext_data:
                        print("   Дополнительно:")
                        for k, v in ext_data.items():
                            print(f"     • {k}: {v}")
                except Exception:
                    pass

                print("\n" + "-" * 50 + "\n")

        except Exception as e:
            print(f"❌ Ошибка поиска: {e}")
            logger.exception("Ошибка в search()")

    # ---------- главное меню ----------
    def print_menu(self):
        print("\n" + "=" * 60)
        print(" PhyloResolve — консольная версия")
        print("=" * 60)
        self.update_status()
        print("""
1) 🌐  Загрузить дерево YFull из сети
2) 📥  Загрузить JSON
3) 📊  Загрузить CSV
4) 🔄  Конвертировать файл
5) 🔎  Поиск
6) 🗑   Очистить базу
0) 🚪  Выход
""")

    def run(self):
        while True:
            self.print_menu()
            choice = input("Выберите пункт меню: ").strip()
            if choice == "1":
                self.load_yfull_online()
            elif choice == "2":
                self.load_file("json")
            elif choice == "3":
                self.load_file("csv")
            elif choice == "4":
                self.convert_file_cli()
            elif choice == "5":
                self.search()
            elif choice == "6":
                self.clear_db()
            elif choice == "0":
                print("До свидания.")
                break
            else:
                print("Неизвестный пункт меню, попробуйте снова.")


if __name__ == "__main__":
    try:
        Cli().run()
    except KeyboardInterrupt:
        print("\n⏹️ Прервано пользователем.")
    except Exception:
        logger.exception("Необработанная ошибка")
        raise
    finally:
        try:
            input("\nНажмите Enter, чтобы закрыть окно...")
        except (EOFError, KeyboardInterrupt):
            pass
