#!/usr/bin/env python3
"""
format_references.py  (standard library only)
Formats a numbered reference list from IDs; bibliographic data always comes from the source APIs.

  python3 scripts/format_references.py PMID:34584091 DOI:10.1038/s41467-021-25960-2 SW:scanpy
  python3 scripts/format_references.py --file ids.txt        # one ID per line ('#' comments allowed)
  options: --refresh (ignore cache), --offline (cache only; error if missing), --no-number

ID kinds
  PMID:<n>  NCBI E-utilities esummary (db=pubmed, retmode=json); <= 3 requests/second
  DOI:<doi> Crossref api.crossref.org/works/<doi>
  ARXIV:<id> arXiv API export.arxiv.org/api/query?id_list=<id> (raw Atom XML is cached; plain-text preprint citation, no URL)
  NCT:<id>  ClinicalTrials.gov API v2 clinicaltrials.gov/api/v2/studies/<id> (selected fields; plain-text registry citation, no URL)
  SW:<key>  entry of review/software_refs.json ("paper": "PMID:..." / "DOI:..." is formatted as the original paper)
Contact: REF_MAILTO=<e-mail> adds a mailto to the User-Agent (Crossref/NCBI etiquette); unset = no address is sent.
Cache: review/references_cache.json (ID -> raw API metadata); a re-run with a cache makes no network calls.
Overrides: review/reference_overrides.json (ID -> corrected fields + "source"/"reason"), applied after the cache read; the cache keeps the raw API response.
  Supported fields (DOI records): "author" (Crossref-style [{"given","family"}, ...]) and "year" (int). Unknown fields are an error.
Output: Markdown numbered list, no hyperlinks (build assertion 5 forbids href="https://...").
Failures: the failing IDs are reported on stderr and the exit code is 1 (hand-written entries are never produced).
"""
import json, os, re, sys, time, copy, urllib.request, urllib.parse, urllib.error
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "review", "references_cache.json")
SWREFS = os.path.join(ROOT, "review", "software_refs.json")
OVERRIDES = os.path.join(ROOT, "review", "reference_overrides.json")
NCT_FIELDS = "NCTId,BriefTitle,OverallStatus,LeadSponsorName,StartDate,CompletionDate,Phase,WhyStopped"


def _user_agent():
    mailto = os.environ.get("REF_MAILTO")
    return "course-reference-formatter/1.0" + (f" (mailto:{mailto})" if mailto else "")


UA = _user_agent()
_last = [0.0]


def _get(url):
    wait = 0.36 - (time.time() - _last[0])  # NCBI: <= 3 requests/second
    if wait > 0:
        time.sleep(wait)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                _last[0] = time.time()
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503) and attempt < 2:
                time.sleep(2 * (attempt + 1)); continue
            raise
    raise RuntimeError("unreachable")


def _get_text(url, min_gap=3.0):
    """arXiv API: raw Atom text; the API terms ask for >= 3 s between requests."""
    wait = min_gap - (time.time() - _last[0])
    if wait > 0:
        time.sleep(wait)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                _last[0] = time.time()
                return r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503) and attempt < 2:
                time.sleep(5 * (attempt + 1)); continue
            raise
    raise RuntimeError("unreachable")


def fetch(kind, key):
    if kind == "ARXIV":
        txt = _get_text("https://export.arxiv.org/api/query?id_list=" + urllib.parse.quote(key))
        if "<entry>" not in txt or "<title>Error</title>" in txt:
            raise ValueError("no arXiv entry")
        return {"atom": txt}
    if kind == "NCT":
        d = _get("https://clinicaltrials.gov/api/v2/studies/" + urllib.parse.quote(key) + "?fields=" + NCT_FIELDS)
        if d.get("protocolSection", {}).get("identificationModule", {}).get("nctId") != key:
            raise ValueError("no ClinicalTrials.gov record")
        return d
    if kind == "PMID":
        d = _get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + urllib.parse.quote(key))
        rec = d.get("result", {}).get(key)
        if not rec or "error" in rec or not rec.get("title"):
            raise ValueError("no esummary record")
        return rec
    d = _get("https://api.crossref.org/works/" + urllib.parse.quote(key, safe="/"))
    if d.get("status") != "ok":
        raise ValueError("crossref status " + str(d.get("status")))
    return d["message"]


def _authors(names):
    out = names[:3]
    return ", ".join(out) + (", et al" if len(names) > 3 else "")


def fmt_pmid(rec, key):
    names = [a["name"] for a in rec.get("authors", []) if a.get("name")]  # includes collective names
    year = (rec.get("pubdate") or rec.get("epubdate") or "").split(" ")[0]
    vol, iss, pg = rec.get("volume", ""), rec.get("issue", ""), rec.get("pages", "")
    loc = year + (";" + vol if vol else "") + ("(" + iss + ")" if iss else "") + (":" + pg if pg else "")
    if not pg and rec.get("elocationid") and not vol:
        loc = year
    doi = next((x["value"] for x in rec.get("articleids", []) if x.get("idtype") == "doi"), "")
    title = rec["title"].rstrip(". ")
    s = f"{_authors(names)}. {title}. *{rec.get('source', '')}*. {loc}. PMID: {key}."
    if doi:
        s += f" doi: {doi}"
    return s.replace("..", ".")


def _initials(given):
    return "".join(p[0] for p in given.replace("-", " ").split() if p)


def fmt_doi(msg, key):
    names = [(a.get("family", "") + " " + _initials(a.get("given", ""))).strip() for a in msg.get("author", []) if a.get("family") or a.get("name")]
    if not names:
        names = ["(no author listed)"]
    date = (msg.get("issued") or msg.get("published-print") or msg.get("published-online") or {}).get("date-parts", [[""]])[0]
    year = str(date[0]) if date and date[0] else ""
    cont = (msg.get("short-container-title") or msg.get("container-title") or [""])[0]
    vol, iss, pg = msg.get("volume", ""), msg.get("issue", ""), msg.get("page", "") or msg.get("article-number", "")
    loc = year + (";" + vol if vol else "") + ("(" + iss + ")" if iss else "") + (":" + pg if pg else "")
    clean = lambda t: re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", t)).strip()  # Crossref titles may hold tags/newlines
    title = clean((msg.get("title") or [""])[0]).rstrip(". ")
    cont = clean(cont).rstrip(".")
    jr = f" *{cont}*." if cont else ""  # books/preprints may have no container title
    return f"{_authors(names)}. {title}.{jr} {loc}. doi: {msg.get('DOI', key)}".replace("..", ".")


def fmt_arxiv(rec, key):
    ns = {"a": "http://www.w3.org/2005/Atom"}
    e = ET.fromstring(rec["atom"]).find("a:entry", ns)
    if e is None:
        raise ValueError("no arXiv entry")
    names = []
    for a in e.findall("a:author", ns):
        parts = " ".join(a.findtext("a:name", "", ns).split()).split(" ")
        names.append((parts[-1] + " " + _initials(" ".join(parts[:-1]))).strip())
    title = " ".join(e.findtext("a:title", "", ns).split()).rstrip(". ")
    year = e.findtext("a:published", "", ns)[:4]
    return f"{_authors(names)}. {title}. *arXiv* preprint. {year}. arXiv: {key}.".replace("..", ".")


def fmt_nct(rec, key):
    p = rec["protocolSection"]
    title = p["identificationModule"]["briefTitle"].rstrip(". ")
    sponsor = p.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name", "(no sponsor listed)").rstrip(".")
    st = p.get("statusModule", {})
    phases = "/".join(x.replace("PHASE", "") for x in p.get("designModule", {}).get("phases", []) if x.startswith("PHASE"))
    status = st.get("overallStatus", "").replace("_", " ").capitalize()
    start = st.get("startDateStruct", {}).get("date", "")[:4]
    end = st.get("completionDateStruct", {}).get("date", "")[:4]
    period = f"{start}-{end}" if start and end else start
    info = "; ".join(x for x in (f"Phase {phases}" if phases else "", status, period) if x)
    return f"{sponsor}. {title}. *ClinicalTrials.gov* registry ({info}). ClinicalTrials.gov: {key}.".replace("..", ".")


def apply_override(kind, key, rec, ov=None):
    """Return rec with the corrected fields of review/reference_overrides.json[kind:key] applied (deep copy; cache untouched)."""
    ov = load_json(OVERRIDES, {}) if ov is None else ov
    o = ov.get(f"{kind}:{key}")
    if not o:
        return rec
    rec = copy.deepcopy(rec)
    for f, v in o.items():
        if f in ("source", "reason"):
            continue
        if kind != "DOI" or f not in ("author", "year"):
            raise ValueError(f"unsupported override field {kind}:{key}.{f}")
        if f == "author":
            rec["author"] = v
        else:
            rec["issued"] = {"date-parts": [[int(v)]]}
    return rec


def load_json(p, default):
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else default


def format_id(raw, get, sw):
    """Format one ID ('PMID:..', 'DOI:..', 'SW:key') using get(kind, key) -> raw API metadata. Raises on failure."""
    kind, _, key = raw.partition(":")
    kind = kind.upper(); key = key.strip()

    def ref(k, v):
        rec = apply_override(k, v, get(k, v))
        fmt = {"PMID": fmt_pmid, "ARXIV": fmt_arxiv, "NCT": fmt_nct}.get(k, fmt_doi)
        return fmt(rec, v)

    if kind in ("PMID", "DOI", "ARXIV", "NCT"):
        return ref(kind, key)
    if kind == "SW":
        e = sw[key]
        ver = (" " + e["version_checked"]) if e.get("version_checked") else ""  # blank = version not verified
        s = f"{e['name']}{ver} ({e['source']})."
        p = e.get("paper")
        if p:
            pk, _, pv = p.partition(":")
            s += " 원 논문: " + ref(pk.upper(), pv)
        return s
    raise ValueError("unknown ID kind")


def main(argv):
    flags = {a for a in argv if a.startswith("--") and a != "--file"}
    args = [a for a in argv if not a.startswith("--")]
    if "--file" in argv:
        fn = argv[argv.index("--file") + 1]
        args = [l.split("#")[0].strip() for l in open(fn, encoding="utf-8") if l.split("#")[0].strip()]
        args = [a for a in args if a != fn]
    cache = load_json(CACHE, {})
    sw = load_json(SWREFS, {})
    failures, lines = [], []

    def get(kind, key):
        ck = f"{kind}:{key}"
        if ck in cache and "--refresh" not in flags:
            return cache[ck]
        if "--offline" in flags:
            raise LookupError("not in cache (offline)")
        cache[ck] = fetch(kind, key)
        return cache[ck]

    for raw in args:
        try:
            lines.append(format_id(raw, get, sw))
        except Exception as ex:  # report, never fabricate
            failures.append((raw, f"{type(ex).__name__}: {ex}"))
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    json.dump(cache, open(CACHE, "w", encoding="utf-8"), ensure_ascii=False, indent=1, sort_keys=True)
    for i, l in enumerate(lines, 1):
        print(l if "--no-number" in flags else f"{i}. {l}")
    if failures:
        print(f"\n[format_references] {len(failures)} ID(s) failed:", file=sys.stderr)
        for r, e in failures:
            print(f"  {r}: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
