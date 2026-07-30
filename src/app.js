const DNA_BASES = new Set(["A", "C", "G", "T", "N"]);
const FEATURE_PALETTE = {
  promoter: "#e8724a",
  cds: "#6555d8",
  gene: "#6555d8",
  terminator: "#d9a63c",
  origin: "#1a9f91",
  rep_origin: "#1a9f91",
  primer_bind: "#e25e9f",
  misc_feature: "#6e8799",
  marker: "#3a86a8",
};

const BROWSER_STORAGE_KEY = "opengene-studio-project-v1";

const ENZYMES = [
  { name: "EcoRI", site: "GAATTC", cut: 1, accent: "#30a78e" },
  { name: "BamHI", site: "GGATCC", cut: 1, accent: "#e8724a" },
  { name: "HindIII", site: "AAGCTT", cut: 1, accent: "#a870e8" },
  { name: "XhoI", site: "CTCGAG", cut: 1, accent: "#e7b546" },
  { name: "NotI", site: "GCGGCCGC", cut: 2, accent: "#5880e6" },
  { name: "NheI", site: "GCTAGC", cut: 1, accent: "#ea5c7b" },
  { name: "BsaI", site: "GGTCTC", cut: 7, accent: "#467b63" },
];

const COMMON_FEATURES_DB = [
  { name: "AmpR", type: "cds", sequence: "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGCAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGG" },
  { name: "KanR", type: "cds", sequence: "ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGACGGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTTCTTCTGA" },
  { name: "mCherry", type: "cds", sequence: "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA" },
  { name: "T7 promoter", type: "promoter", sequence: "TAATACGACTCACTATAGGG" },
  { name: "T7 terminator", type: "terminator", sequence: "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG" },
  { name: "lac promoter", type: "promoter", sequence: "GCATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGA" },
  { name: "6xHis", type: "cds", sequence: "CATCACCATCACCATCAC" },
  { name: "GFP", type: "cds", sequence: "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA" },
  { name: "pUC ori", type: "origin", sequence: "CCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAG" }
];

const GENETIC_CODE = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", ATG: "M", TAA: "*", TAG: "*", TGA: "*",
  GCT: "A", GCC: "A", GCA: "A", GCG: "A", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
  CCT: "P", CCC: "P", CCA: "P", CCG: "P", AAT: "N", AAC: "N", GAT: "D", GAC: "D",
  CAA: "Q", CAG: "Q", AAA: "K", AAG: "K", TAT: "Y", TAC: "Y", CAT: "H", CAC: "H",
  TGT: "C", TGC: "C", TGG: "W", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
  AGA: "R", AGG: "R", ACT: "T", ACC: "T", ACA: "T", ACG: "T", GTT: "V",
  GTC: "V", GTA: "V", GTG: "V", TCT: "S", TCC: "S", TCA: "S", TCG: "S",
  AGT: "S", AGC: "S", CTT: "L", CTC: "L", CTA: "L", CTG: "L",
  ATT: "I", ATC: "I", ATA: "I", GAA: "E", GAG: "E",
};

const STOP_CODONS = new Set(["TAA", "TAG", "TGA"]);
const PREFERRED_CODON = Object.entries(GENETIC_CODE).reduce((table, [codon, aminoAcid]) => {
  if (!table[aminoAcid]) table[aminoAcid] = codon;
  return table;
}, {});

const state = {
  documents: [],
  activeId: null,
  view: "map",
  mode: "circular",
  selection: { start: 426, end: 426 },
  activeFeatureId: null,
  showSites: true,
  activeTool: "restriction",
  undoStack: [],
  redoStack: [],
  history: [],
  sidebarOpen: false,
  editOperation: "replace",
  darkMap: false,
  search: "",
  reads: [],
  verificationRun: null,
};

function uid(prefix = "id") {
  return `${prefix}-${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
}

function escapeHtml(value = "") {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function normalizeSequence(value) {
  const sequence = String(value || "")
    .toUpperCase()
    .replace(/\s+/g, "")
    .replace(/U/g, "T");
  return [...sequence].filter((base) => DNA_BASES.has(base)).join("");
}

function reverseComplement(sequence) {
  const bases = { A: "T", T: "A", C: "G", G: "C", N: "N" };
  return [...sequence.toUpperCase()].reverse().map((base) => bases[base] || "N").join("");
}

function gcPercent(sequence) {
  if (!sequence.length) return 0;
  return Math.round((([...sequence].filter((base) => base === "G" || base === "C").length / sequence.length) * 1000)) / 10;
}

function baseSequence(length) {
  const alphabet = ["G", "C", "T", "A", "A", "G", "C", "T", "G", "A", "T", "C"];
  let value = "";
  let seed = 937;
  for (let index = 0; index < length; index += 1) {
    seed = (seed * 48271) % 2147483647;
    value += alphabet[(seed + index * 7) % alphabet.length];
  }
  const motifs = [
    [188, "GAATTC"], [1042, "GGATCC"], [1851, "GCTAGC"], [2472, "CTCGAG"],
    [3287, "AAGCTT"], [4159, "GGTCTC"], [5016, "GCGGCCGC"],
  ];
  for (const [offset, motif] of motifs) {
    value = `${value.slice(0, offset)}${motif}${value.slice(offset + motif.length)}`;
  }
  return value;
}

function createDemoDocument() {
  return {
    id: uid("doc"),
    name: "pOGS–Aurora",
    description: "Starter construct · local-only demonstration project",
    sequence: baseSequence(5821),
    topology: "circular",
    moleculeType: "DNA",
    revision: 1,
    metadata: { source: "OpenGene Starter Library", createdAt: new Date().toISOString() },
    features: [
      feature("promoter", "T7 promoter", 146, 202, 1, { note: "Inducible transcription start" }),
      feature("primer_bind", "OGS-Fwd", 364, 387, 1, { label: "OGS-Fwd" }),
      feature("cds", "mNeonGreen", 521, 1244, 1, { gene: "mNeonGreen", translation: "fluorescent reporter" }),
      feature("terminator", "rrnB T1 terminator", 1279, 1374, 1, {}),
      feature("marker", "Kanamycin resistance", 1714, 2525, -1, { gene: "kanR" }),
      feature("origin", "pUC origin", 3273, 4095, 1, { note: "High-copy replication origin" }),
      feature("primer_bind", "OGS-Rev", 4978, 5003, -1, { label: "OGS-Rev" }),
      feature("misc_feature", "Assembly junction", 5734, 85, 1, { note: "Origin-spanning linker", wrapsOrigin: true }),
    ],
  };
}

function mutateBase(sequence, position) {
  const bases = ["A", "C", "G", "T"];
  const existing = sequence[position];
  const replacement = bases.find((base) => base !== existing) || "N";
  return `${sequence.slice(0, position)}${replacement}${sequence.slice(position + 1)}`;
}

function createDemoReads(document) {
  const makeRead = (name, start, length, mutationOffset = null) => {
    let sequence = document.sequence.slice(start, start + length);
    if (mutationOffset !== null && mutationOffset < sequence.length) sequence = mutateBase(sequence, mutationOffset);
    return {
      id: uid("read"),
      name,
      sequence,
      source: "starter read",
      quality: Array.from({ length: sequence.length }, (_, index) => Math.max(22, 39 - (index % 11))),
      createdAt: new Date().toISOString(),
    };
  };
  return [
    makeRead("Aurora forward read", 300, 170),
    makeRead("Aurora reporter read", 580, 180, 74),
    makeRead("Aurora backbone read", 3250, 160),
  ];
}

function feature(type, name, start, end, strand = 1, qualifiers = {}) {
  return {
    id: uid("feature"),
    type,
    name,
    start: Number(start),
    end: Number(end),
    strand: Number(strand) >= 0 ? 1 : -1,
    wrapsOrigin: Boolean(qualifiers.wrapsOrigin || start > end),
    qualifiers,
  };
}

function activeDocument() {
  return state.documents.find((document) => document.id === state.activeId) || state.documents[0];
}

function clone(value) {
  return JSON.parse(JSON.stringify(value));
}

function selectionLength(document = activeDocument()) {
  const { start, end } = state.selection;
  if (!document || start === end) return 0;
  return end >= start ? end - start : document.sequence.length - start + end;
}

function selectionSequence(document = activeDocument()) {
  if (!document) return "";
  const { start, end } = state.selection;
  if (start === end) return "";
  if (end > start) return document.sequence.slice(start, end);
  return `${document.sequence.slice(start)}${document.sequence.slice(0, end)}`;
}

function featureSegments(item, sequenceLength) {
  if (item.wrapsOrigin || item.start > item.end) {
    return [[item.start, sequenceLength], [0, item.end]];
  }
  return [[item.start, item.end]];
}

function featureLength(item, document = activeDocument()) {
  return featureSegments(item, document.sequence.length).reduce((sum, [start, end]) => sum + Math.max(0, end - start), 0);
}

function categoryColor(item) {
  return FEATURE_PALETTE[item.type] || FEATURE_PALETTE.misc_feature;
}

function formatBp(value) {
  return `${new Intl.NumberFormat().format(value)} bp`;
}

function formatPosition(position) {
  return new Intl.NumberFormat().format(Math.max(0, position) + 1);
}

function plural(value, unit) {
  return `${value} ${unit}${value === 1 ? "" : "s"}`;
}

function locationLabel(item, document = activeDocument()) {
  const range = item.wrapsOrigin || item.start > item.end
    ? `${formatPosition(item.start)} → ${formatPosition(item.end)} (wrap)`
    : `${formatPosition(item.start)} – ${formatPosition(Math.max(item.start, item.end - 1))}`;
  return `${range} · ${item.strand === 1 ? "+" : "−"}`;
}

function makeSnapshot() {
  return {
    documents: clone(state.documents),
    activeId: state.activeId,
    selection: clone(state.selection),
    activeFeatureId: state.activeFeatureId,
  };
}

function restoreSnapshot(snapshot) {
  state.documents = clone(snapshot.documents);
  state.activeId = snapshot.activeId;
  state.selection = clone(snapshot.selection);
  state.activeFeatureId = snapshot.activeFeatureId;
}

function commit(description, mutation, detail = "") {
  const before = makeSnapshot();
  mutation();
  const document = activeDocument();
  if (document) document.revision = (document.revision || 0) + 1;
  state.undoStack.push({ description, before, after: makeSnapshot() });
  state.redoStack = [];
  recordHistory({ type: "edit", title: description, detail: detail || `${document?.name || "Document"} · revision ${document?.revision || 1}` });
  render();
  showToast(description);
}

function undo() {
  const command = state.undoStack.pop();
  if (!command) return showToast("Nothing to undo");
  restoreSnapshot(command.before);
  state.redoStack.push(command);
  persistLocalProject();
  render();
  showToast(`Undid: ${command.description}`);
}

function redo() {
  const command = state.redoStack.pop();
  if (!command) return showToast("Nothing to redo");
  restoreSnapshot(command.after);
  state.undoStack.push(command);
  persistLocalProject();
  render();
  showToast(`Redid: ${command.description}`);
}

function remapAfterEdit(document, position, deletedLength, insertedLength) {
  const delta = insertedLength - deletedLength;
  const deletionEnd = position + deletedLength;
  const mapPoint = (value, affinity) => {
    if (value < position) return value;
    if (value > deletionEnd) return value + delta;
    if (value === position && affinity === "start") return position;
    return position + insertedLength;
  };
  for (const item of document.features) {
    if (item.wrapsOrigin) {
      item.start = Math.max(0, Math.min(document.sequence.length, mapPoint(item.start, "start")));
      item.end = Math.max(0, Math.min(document.sequence.length, mapPoint(item.end, "end")));
      continue;
    }
    const oldStart = item.start;
    const oldEnd = item.end;
    item.start = Math.max(0, mapPoint(oldStart, "start"));
    item.end = Math.max(item.start, mapPoint(oldEnd, "end"));
    if (oldStart < position && oldEnd > deletionEnd) item.end += 0;
  }
}

function replaceSequenceRegion(start, end, replacement) {
  const document = activeDocument();
  if (!document) return;
  const value = normalizeSequence(replacement);
  if (start < 0 || end < start || end > document.sequence.length) return;
  const previous = document.sequence.slice(start, end);
  commit(
    previous.length === 0 ? `Inserted ${formatBp(value.length)}` : value.length === 0 ? `Deleted ${formatBp(previous.length)}` : "Replaced sequence region",
    () => {
      document.sequence = `${document.sequence.slice(0, start)}${value}${document.sequence.slice(end)}`;
      remapAfterEdit(document, start, previous.length, value.length);
      state.selection = { start, end: start + value.length };
    },
    `Edited bases ${formatPosition(start)}–${formatPosition(Math.max(start, end - 1))}`,
  );
}

function addFeatureFromForm(values) {
  const document = activeDocument();
  const start = Math.max(0, Number(values.start) - 1);
  const end = Math.max(0, Number(values.end));
  if (!values.name.trim() || !Number.isFinite(start) || !Number.isFinite(end)) return showToast("Enter a name and valid coordinates");
  if (start > document.sequence.length || end > document.sequence.length || start === end) return showToast("Feature coordinates are outside this sequence");
  const item = feature(values.type, values.name.trim(), start, end, values.strand, { note: values.note.trim(), wrapsOrigin: start > end });
  commit(`Added ${item.name}`, () => {
    document.features.push(item);
    state.activeFeatureId = item.id;
    state.selection = { start: item.start, end: item.end };
  }, `Annotated ${formatBp(featureLength(item, document))}`);
  closeModal();
}

function deleteFeature(id) {
  const document = activeDocument();
  const item = document.features.find((candidate) => candidate.id === id);
  if (!item) return;
  commit(`Removed ${item.name}`, () => {
    document.features = document.features.filter((candidate) => candidate.id !== id);
    state.activeFeatureId = null;
  });
}

function selectFeature(id) {
  const document = activeDocument();
  const item = document.features.find((candidate) => candidate.id === id);
  if (!item) return;
  state.activeFeatureId = id;
  state.selection = { start: item.start, end: item.end };
  render();
}

function openDocument(id) {
  const target = state.documents.find((document) => document.id === id);
  if (!target) return showToast("That document is no longer available");
  state.activeId = target.id;
  state.selection = { start: 0, end: 0 };
  state.activeFeatureId = null;
  state.view = "map";
  persistLocalProject();
  render();
}

function duplicateActiveDocument() {
  const source = activeDocument();
  const copy = clone(source);
  copy.id = uid("doc");
  copy.name = `${source.name} copy`;
  copy.description = `Local duplicate of ${source.name}`;
  copy.revision = 1;
  copy.metadata = { ...copy.metadata, parentDocumentId: source.id, duplicated: true };
  state.documents.push(copy);
  state.activeId = copy.id;
  recordOperation({ type: "edit", title: "Duplicated document", detail: `${source.name} → ${copy.name}`, inputs: [source.id], output: copy.id });
  state.view = "map";
  render();
  showToast(`Created ${copy.name}`);
}

function findSites(sequence, enzyme, isCircular = false) {
  const hits = [];
  const searchable = isCircular ? sequence + sequence.slice(0, Math.max(0, enzyme.site.length - 1)) : sequence;
  const forwardSite = enzyme.site;
  const reverseSite = reverseComplement(enzyme.site);

  let index = searchable.indexOf(forwardSite);
  while (index !== -1) {
    const pos = index % sequence.length;
    if (!hits.some(h => h.position === pos)) {
      hits.push({ ...enzyme, position: pos, cutPosition: (index + enzyme.cut) % sequence.length });
    }
    index = searchable.indexOf(forwardSite, index + 1);
  }

  if (forwardSite !== reverseSite) {
    let revIndex = searchable.indexOf(reverseSite);
    while (revIndex !== -1) {
      const pos = revIndex % sequence.length;
      if (!hits.some(h => h.position === pos)) {
        const cutOffset = enzyme.site.length - enzyme.cut;
        hits.push({ ...enzyme, position: pos, cutPosition: (revIndex + cutOffset + sequence.length) % sequence.length });
      }
      revIndex = searchable.indexOf(reverseSite, revIndex + 1);
    }
  }

  return hits;
}

function allSites(document = activeDocument()) {
  if (!document) return [];
  return ENZYMES.flatMap((enzyme) => findSites(document.sequence, enzyme, document.topology === "circular"));
}

function digestFragments(document, siteList) {
  const cuts = [...new Set(siteList.map((site) => site.cutPosition).filter((position) => position >= 0 && position < document.sequence.length))].sort((a, b) => a - b);
  if (!cuts.length) return [{ start: 0, end: document.sequence.length, length: document.sequence.length }];
  if (document.topology === "linear") {
    const points = [0, ...cuts, document.sequence.length];
    return points.slice(0, -1).map((start, index) => ({ start, end: points[index + 1], length: points[index + 1] - start }));
  }
  if (cuts.length === 1) return [{ start: cuts[0], end: cuts[0], length: document.sequence.length }];
  return cuts.map((start, index) => {
    const end = cuts[(index + 1) % cuts.length];
    return { start, end, length: end > start ? end - start : document.sequence.length - start + end };
  });
}

const NN_PARAMS = {
  AA: { H: -7.9, S: -22.2 }, TT: { H: -7.9, S: -22.2 },
  AT: { H: -7.2, S: -20.4 }, TA: { H: -7.2, S: -21.3 },
  CA: { H: -8.5, S: -22.7 }, TG: { H: -8.5, S: -22.7 },
  GT: { H: -8.4, S: -22.4 }, AC: { H: -8.4, S: -22.4 },
  CT: { H: -7.8, S: -21.0 }, AG: { H: -7.8, S: -21.0 },
  GA: { H: -8.2, S: -22.2 }, TC: { H: -8.2, S: -22.2 },
  CG: { H: -10.6, S: -27.2 }, GC: { H: -9.8, S: -24.4 },
  GG: { H: -8.0, S: -19.9 }, CC: { H: -8.0, S: -19.9 },
};

function calculateTm(sequence, naConc = 0.05, primerConc = 5e-7) {
  const normalized = normalizeSequence(sequence);
  if (normalized.length < 14) {
    const at = [...normalized].filter((base) => base === "A" || base === "T").length;
    const gc = [...normalized].filter((base) => base === "G" || base === "C").length;
    return at * 2 + gc * 4;
  }
  let dH = 0;
  let dS = 0;
  let validPairs = 0;
  for (let i = 0; i < normalized.length - 1; i++) {
    const pair = normalized.slice(i, i + 2);
    if (NN_PARAMS[pair]) {
      dH += NN_PARAMS[pair].H;
      dS += NN_PARAMS[pair].S;
      validPairs++;
    }
  }
  if (validPairs === 0) return 0;
  const first = normalized[0];
  const last = normalized[normalized.length - 1];
  if (first === "G" || first === "C") { dH += 0.1; dS += -2.8; }
  else if (first === "A" || first === "T") { dH += 2.3; dS += 4.1; }
  if (last === "G" || last === "C") { dH += 0.1; dS += -2.8; }
  else if (last === "A" || last === "T") { dH += 2.3; dS += 4.1; }

  dS += 0.368 * (normalized.length - 1) * Math.log(naConc);
  const tm = (dH * 1000) / (dS + 1.987 * Math.log(primerConc / 4)) - 273.15;
  return Math.max(0, Math.round(tm * 10) / 10);
}

function primerCandidates(document = activeDocument()) {
  if (!document) return [];
  const { start, end } = state.selection;
  const anchorStart = Math.max(0, Math.min(document.sequence.length - 22, start || 360));
  const anchorEnd = Math.max(22, Math.min(document.sequence.length, (end || 5003)));
  const forward = document.sequence.slice(anchorStart, anchorStart + 22);
  const reverse = reverseComplement(document.sequence.slice(anchorEnd - 22, anchorEnd));
  const product = end > start ? end - start : Math.max(120, anchorEnd - anchorStart);
  return [
    { name: "OGS-forward", direction: "Forward", sequence: forward, position: anchorStart, tm: calculateTm(forward), gc: gcPercent(forward), product },
    { name: "OGS-reverse", direction: "Reverse", sequence: reverse, position: Math.max(0, anchorEnd - 22), tm: calculateTm(reverse), gc: gcPercent(reverse), product },
  ];
}

function bestUngappedAlignment(reference, inputSequence, circular = false) {
  const sequence = normalizeSequence(inputSequence);
  if (!reference?.length || !sequence.length) return null;
  const orientations = [
    { sequence, strand: "+" },
    { sequence: reverseComplement(sequence), strand: "−" },
  ];
  const searchable = circular ? `${reference}${reference.slice(0, Math.max(0, sequence.length - 1))}` : reference;
  let best = null;
  for (const orientation of orientations) {
    const maximum = Math.max(0, searchable.length - orientation.sequence.length);
    for (let offset = 0; offset <= maximum; offset += 1) {
      if (!circular && offset + orientation.sequence.length > reference.length) continue;
      let matches = 0;
      const variants = [];
      for (let index = 0; index < orientation.sequence.length; index += 1) {
        const position = (offset + index) % reference.length;
        const base = orientation.sequence[index];
        const referenceBase = reference[position];
        if (base === referenceBase || base === "N") matches += 1;
        else variants.push({ position, reference: referenceBase, alternate: base });
      }
      const candidate = { start: offset % reference.length, end: (offset + orientation.sequence.length) % reference.length, strand: orientation.strand, sequence: orientation.sequence, matches, variants, identity: matches / orientation.sequence.length };
      if (!best || candidate.matches > best.matches || (candidate.matches === best.matches && candidate.variants.length < best.variants.length)) best = candidate;
    }
  }
  return best;
}

function verificationSummary(document = activeDocument()) {
  if (!document) return { alignments: [], coverage: [], variants: [], meanIdentity: 0, consensus: "" };
  const coverage = Array.from({ length: document.sequence.length }, () => ({ depth: 0, bases: {} }));
  const alignments = state.reads.map((read) => {
    const alignment = bestUngappedAlignment(document.sequence, read.sequence, document.topology === "circular");
    if (!alignment) return { read, alignment: null };
    for (let index = 0; index < alignment.sequence.length; index += 1) {
      const position = (alignment.start + index) % document.sequence.length;
      const base = alignment.sequence[index];
      if (base === "N") continue;
      coverage[position].depth += 1;
      coverage[position].bases[base] = (coverage[position].bases[base] || 0) + 1;
    }
    return { read, alignment };
  });
  const variants = [];
  const consensus = [...document.sequence];
  for (let position = 0; position < coverage.length; position += 1) {
    const entry = coverage[position];
    const calls = Object.entries(entry.bases).sort(([, a], [, b]) => b - a);
    if (!calls.length) continue;
    const [base, count] = calls[0];
    const confident = count / entry.depth >= .7;
    if (confident) consensus[position] = base;
    if (confident && base !== document.sequence[position]) variants.push({ position, reference: document.sequence[position], alternate: base, depth: entry.depth, support: count });
  }
  const mapped = alignments.filter((item) => item.alignment);
  return {
    alignments,
    coverage,
    variants,
    consensus: consensus.join(""),
    meanIdentity: mapped.length ? mapped.reduce((sum, item) => sum + item.alignment.identity, 0) / mapped.length : 0,
    coveredBases: coverage.filter((item) => item.depth > 0).length,
  };
}

function coverageBins(summary, binCount = 72) {
  if (!summary.coverage?.length) return [];
  const binSize = Math.ceil(summary.coverage.length / binCount);
  return Array.from({ length: Math.ceil(summary.coverage.length / binSize) }, (_, bin) => {
    const items = summary.coverage.slice(bin * binSize, (bin + 1) * binSize);
    return items.reduce((total, item) => total + item.depth, 0) / Math.max(1, items.length);
  });
}

function confidenceLabel(identity) {
  if (identity >= .995) return "validated candidate";
  if (identity >= .97) return "review candidate";
  return "needs review";
}

function recordHistory(entry) {
  state.history.unshift({
    id: uid("event"),
    timestamp: new Date().toISOString(),
    inputs: [activeDocument()?.id].filter(Boolean),
    output: activeDocument()?.id || null,
    warnings: [],
    ...entry,
  });
  state.history = state.history.slice(0, 24);
  persistLocalProject();
}

function createInitialHistory() {
  state.history = [
    {
      id: uid("event"), type: "import", title: "Created starter construct", detail: "pOGS–Aurora · 5,821 bp circular DNA", timestamp: new Date(Date.now() - 1000 * 60 * 48).toISOString(), inputs: [], output: state.activeId, warnings: [],
    },
    {
      id: uid("event"), type: "annotation", title: "Applied starter annotations", detail: "8 features from OpenGene starter library", timestamp: new Date(Date.now() - 1000 * 60 * 47).toISOString(), inputs: [state.activeId], output: state.activeId, warnings: [],
    },
  ];
}

function parseFasta(text, filename = "Imported sequence") {
  const lines = text.replace(/\r/g, "").split("\n");
  const header = lines.find((line) => line.startsWith(">")) || `>${filename.replace(/\.[^.]+$/, "")}`;
  const name = header.slice(1).trim().split(/\s{2,}|\|/)[0] || filename;
  const description = header.slice(1).trim();
  const sequence = normalizeSequence(lines.filter((line) => !line.startsWith(">")).join(""));
  if (!sequence) throw new Error("No DNA sequence found in this FASTA file.");
  return { id: uid("doc"), name, description, sequence, topology: "linear", moleculeType: "DNA", revision: 1, metadata: { source: filename, importedAt: new Date().toISOString() }, features: [] };
}

function parseGenBank(text, filename = "Imported sequence") {
  const locus = text.match(/^LOCUS\s+([^\s]+)\s+(\d+)\s+bp.*?(circular|linear)?/im);
  const name = locus?.[1] || filename.replace(/\.[^.]+$/, "");
  const topology = locus?.[3]?.toLowerCase() || "linear";
  const origin = text.search(/^ORIGIN\b/im);
  if (origin < 0) throw new Error("This GenBank file does not include an ORIGIN sequence section.");
  const sequence = normalizeSequence(text.slice(origin).replace(/^.*?\n/, "").replace(/\/\//g, ""));
  if (!sequence) throw new Error("No DNA sequence found in this GenBank file.");
  const features = [];
  const lines = text.replace(/\r/g, "").split("\n");
  let current = null;
  for (const line of lines) {
    const locationMatch = line.match(/^\s{5}(\S+)\s+(.+)$/);
    if (locationMatch && locationMatch[1] !== "source") {
      if (current) features.push(current);
      const [, type, rawLocation] = locationMatch;
      const complement = /complement/i.test(rawLocation);
      const numbers = [...rawLocation.matchAll(/(\d+)/g)].map((match) => Number(match[1]));
      if (numbers.length >= 2) {
        const start = Math.max(0, numbers[0] - 1);
        const end = numbers[1];
        current = feature(type.toLowerCase(), type, start, end, complement ? -1 : 1, {});
      } else current = null;
      continue;
    }
    if (current) {
      const qualifier = line.match(/^\s+\/(label|gene|note)="?([^"\n]+)"?/);
      if (qualifier) {
        const [, key, value] = qualifier;
        current.qualifiers[key] = value.trim();
        if (["label", "gene"].includes(key)) current.name = value.trim();
      }
    }
  }
  if (current) features.push(current);
  return { id: uid("doc"), name, description: `Imported from ${filename}`, sequence, topology, moleculeType: "DNA", revision: 1, metadata: { source: filename, importedAt: new Date().toISOString() }, features };
}

function parseFastq(text, filename = "Imported reads") {
  const lines = text.replace(/\r/g, "").split("\n").filter((line) => line.length > 0);
  if (!lines.length || lines.length % 4 !== 0) throw new Error("FASTQ records must contain complete four-line entries.");
  const reads = [];
  for (let index = 0; index < lines.length; index += 4) {
    const [header, rawSequence, separator, rawQuality] = lines.slice(index, index + 4);
    if (!header.startsWith("@") || !separator.startsWith("+")) throw new Error(`Invalid FASTQ record at line ${index + 1}.`);
    const sequence = normalizeSequence(rawSequence);
    if (!sequence || rawQuality.length !== rawSequence.trim().length) throw new Error(`FASTQ sequence/quality length mismatch at line ${index + 1}.`);
    reads.push({ id: uid("read"), name: header.slice(1).trim() || `${filename} read ${reads.length + 1}`, sequence, source: filename, quality: [...rawQuality].map((char) => Math.max(0, char.charCodeAt(0) - 33)), createdAt: new Date().toISOString() });
  }
  return reads;
}

function fastaText(document) {
  const header = `>${document.name}${document.description ? ` ${document.description}` : ""}`;
  const chunks = document.sequence.match(/.{1,70}/g) || [];
  return `${header}\n${chunks.join("\n")}\n`;
}

function genbankLocation(item) {
  const range = item.wrapsOrigin || item.start > item.end
    ? `join(${item.start + 1}..${activeDocument().sequence.length},1..${item.end})`
    : `${item.start + 1}..${item.end}`;
  return item.strand === -1 ? `complement(${range})` : range;
}

function genbankText(document) {
  const date = new Date().toISOString().slice(0, 10).replaceAll("-", "-");
  const lines = [
    `LOCUS       ${document.name.slice(0, 16).padEnd(16)} ${String(document.sequence.length).padStart(7)} bp    DNA     ${document.topology.padEnd(8)} ${date}`,
    `DEFINITION  ${document.description || "OpenGene Studio sequence"}.`,
    "FEATURES             Location/Qualifiers",
    `     source          1..${document.sequence.length}`,
    `                     /organism="synthetic DNA construct"`,
  ];
  for (const item of document.features) {
    lines.push(`     ${item.type.slice(0, 15).padEnd(15)} ${genbankLocation(item)}`);
    lines.push(`                     /label="${String(item.name).replaceAll('"', "'")}"`);
    for (const [key, value] of Object.entries(item.qualifiers || {})) {
      if (key === "wrapsOrigin") continue;
      lines.push(`                     /${key}="${String(value).replaceAll('"', "'")}"`);
    }
  }
  lines.push("ORIGIN");
  for (let offset = 0; offset < document.sequence.length; offset += 60) {
    const block = document.sequence.slice(offset, offset + 60).toLowerCase().match(/.{1,10}/g)?.join(" ") || "";
    lines.push(`${String(offset + 1).padStart(9)} ${block}`);
  }
  lines.push("//", "");
  return lines.join("\n");
}

function download(filename, data, type = "text/plain") {
  const blob = new Blob([data], { type });
  const href = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = href;
  anchor.download = filename;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  setTimeout(() => URL.revokeObjectURL(href), 2000);
}

function projectSnapshot() {
  return {
    format: "OpenGene Studio browser project",
    formatVersion: 1,
    appVersion: "0.1.0",
    savedAt: new Date().toISOString(),
    activeId: state.activeId,
    documents: state.documents,
    history: state.history,
    reads: state.reads,
    verificationRun: state.verificationRun,
  };
}

function persistLocalProject() {
  try {
    if (!state.documents.length) return;
    localStorage.setItem(BROWSER_STORAGE_KEY, JSON.stringify(projectSnapshot()));
  } catch {
  }
}

function restoreLocalProject() {
  try {
    const raw = localStorage.getItem(BROWSER_STORAGE_KEY);
    if (!raw) return false;
    const project = JSON.parse(raw);
    if (!Array.isArray(project.documents) || !project.documents.length) return false;
    state.documents = project.documents;
    state.activeId = project.activeId || project.documents[0].id;
    state.history = Array.isArray(project.history) ? project.history : [];
    state.reads = Array.isArray(project.reads) ? project.reads : [];
    state.verificationRun = project.verificationRun || null;
    state.selection = { start: 0, end: 0 };
    return true;
  } catch {
    return false;
  }
}

function render() {
  const host = runtimeHost();
  host.innerHTML = renderWorkspace();
  updateShell();
  requestAnimationFrame(() => {
    attachInteractiveMap();
    attachSequenceTouchSelection();
  });
}

function runtimeHost() {
  const workspace = document.querySelector("#workspace") || document.querySelector("#app") || document.body;
  let host = workspace.querySelector?.("#runtime-content");
  if (host) return host;
  host = document.createElement("div");
  host.id = "runtime-content";
  host.className = "runtime-content";
  const grid = workspace.querySelector?.(".workspace-grid");
  if (grid) {
    grid.replaceChildren(host);
  } else {
    workspace.append(host);
  }
  return host;
}

function renderWorkspace() {
  const document = activeDocument();
  if (!document) return renderEmptyState();
  if (state.view === "projects") {
    return `${renderProjectsView(document)}<div id="modal-layer" aria-live="polite"></div><div id="toast" role="status" aria-live="polite"></div>`;
  }
  const title = `${escapeHtml(document.name)} <span class="revision-badge">r${document.revision || 1}</span>`;
  const description = escapeHtml(document.description || "Local sequence document");
  return `
    <section class="workspace-header">
      <div class="workspace-heading">
        <div class="eyebrow"><span class="status-dot"></span> LOCAL PROJECT · ${document.topology.toUpperCase()}</div>
        <h1>${title}</h1>
        <p>${description}</p>
      </div>
      <div class="workspace-actions" aria-label="Document actions">
        <button class="icon-action" data-action="toggle-sites" aria-pressed="${state.showSites}" title="Toggle restriction sites">${icon("scissors")}</button>
        <button class="soft-action" data-action="edit-sequence">${icon("edit")} <span>Edit bases</span></button>
        <button class="primary-action" data-action="add-feature">${icon("plus")} <span>Add annotation</span></button>
      </div>
    </section>
    ${renderSelectionBar(document)}
    <section class="workspace-content view-${state.view}">
      ${state.view === "map" ? renderMapView(document) : ""}
      ${state.view === "sequence" ? renderSequenceView(document) : ""}
      ${state.view === "features" ? renderFeaturesView(document) : ""}
      ${state.view === "tools" ? renderToolsView(document) : ""}
      ${state.view === "history" ? renderHistoryView(document) : ""}
    </section>
    <div id="modal-layer" aria-live="polite"></div>
    <div id="toast" role="status" aria-live="polite"></div>
  `;
}

function renderProjectsView(document) {
  const totalBases = state.documents.reduce((sum, item) => sum + item.sequence.length, 0);
  const candidateCount = state.documents.filter((item) => item.metadata?.candidate || item.metadata?.consensus).length;
  return `<section class="workspace-header projects-header"><div class="workspace-heading"><div class="eyebrow"><span class="status-dot"></span> LOCAL-FIRST PROJECT SPACE</div><h1>Project library</h1><p>${state.documents.length} document${state.documents.length === 1 ? "" : "s"} · ${new Intl.NumberFormat().format(totalBases)} total bp · no cloud account required</p></div><div class="workspace-actions"><button class="soft-action" data-action="import">${icon("upload")} Import files</button><button class="primary-action" data-action="save-project">${icon("download")} Save project</button></div></section><section class="project-summary-grid"><article class="project-summary-card"><span class="section-label">DOCUMENTS</span><strong>${state.documents.length}</strong><small>plasmids, fragments, candidates</small></article><article class="project-summary-card"><span class="section-label">PROVENANCE</span><strong>${state.history.length}</strong><small>immutable local operation records</small></article><article class="project-summary-card"><span class="section-label">VERIFICATION</span><strong>${state.reads.length}</strong><small>sequence reads available for alignment</small></article><article class="project-summary-card"><span class="section-label">CANDIDATES</span><strong>${candidateCount}</strong><small>assembly or consensus output documents</small></article></section><section class="surface project-library-surface"><div class="surface-heading"><div><span class="section-label">DOCUMENT COLLECTION</span><h2>Open a construct</h2></div><button class="soft-action" data-action="duplicate-document">${icon("copy")} Duplicate active</button></div><div class="project-card-grid">${state.documents.map((item) => `<button class="project-document-card ${item.id === document.id ? "is-active" : ""}" data-action="open-document" data-document-id="${item.id}"><span class="project-document-icon">${icon(item.metadata?.consensus ? "verify" : item.metadata?.candidate ? "assembly" : "helix")}</span><span class="project-document-copy"><strong>${escapeHtml(item.name)}</strong><small>${new Intl.NumberFormat().format(item.sequence.length)} bp · ${escapeHtml(item.topology)} · r${item.revision || 1}</small><em>${item.metadata?.consensus ? "consensus candidate" : item.metadata?.candidate ? "assembly candidate" : item.id === document.id ? "active" : "sequence document"}</em></span><span class="project-card-arrow">${icon("arrow")}</span></button>`).join("")}</div></section><section class="project-integrity"><div><span class="project-integrity-icon">${icon("check")}</span><div><strong>Portable project snapshot ready</strong><small>Manifest, documents, provenance, and locally added read sequences are included in the downloaded <code>.ogs.json</code> snapshot.</small></div></div><button class="text-button" data-action="export-manifest">Inspect manifest</button></section>`;
}

function renderEmptyState() {
  return `<section class="empty-state"><div class="empty-orb">${icon("helix")}</div><h1>Start a local design</h1><p>Import a FASTA or GenBank sequence to begin.</p><button class="primary-action" data-action="import">${icon("upload")} Import sequence</button></section>`;
}

function renderSelectionBar(document) {
  const size = selectionLength(document);
  const selected = document.features.find((item) => item.id === state.activeFeatureId);
  const label = selected
    ? `${escapeHtml(selected.name)} · ${locationLabel(selected, document)}`
    : size ? `${formatPosition(state.selection.start)}–${formatPosition((state.selection.end || 1) - 1)} · ${formatBp(size)}` : `Cursor at ${formatPosition(state.selection.start)}`;
  return `<div class="selection-bar">
    <span class="selection-icon">${selected ? icon("tag") : icon("cursor")}</span>
    <span class="selection-label">${label}</span>
    <span class="selection-grow"></span>
    ${size ? `<button class="compact-action" data-action="copy-selection">Copy ${formatBp(size)}</button>` : ""}
    <button class="compact-action" data-action="zoom-selection">${icon("focus")} Focus</button>
  </div>`;
}

function renderMapView(document) {
  const sites = allSites(document);
  const visibleFeatures = document.features.slice().sort((a, b) => featureLength(b, document) - featureLength(a, document));
  return `
    <div class="map-layout">
      <section class="surface map-surface">
        <div class="surface-heading map-heading">
          <div>
            <span class="section-label">DESIGN CANVAS</span>
            <h2>${state.mode === "circular" ? "Circular construct map" : "Linear construct map"}</h2>
          </div>
          <div class="segmented-control" role="group" aria-label="Map mode">
            <button data-action="set-map-mode" data-mode="circular" class="${state.mode === "circular" ? "is-active" : ""}">${icon("circle")} Circular</button>
            <button data-action="set-map-mode" data-mode="linear" class="${state.mode === "linear" ? "is-active" : ""}">${icon("line")} Linear</button>
          </div>
        </div>
        <div class="map-canvas ${state.darkMap ? "dark-map" : ""}" id="map-canvas" aria-label="Interactive plasmid map">
          ${state.mode === "circular" ? renderCircularMap(document, visibleFeatures, sites) : renderLinearMap(document, visibleFeatures, sites)}
        </div>
        <div class="map-controls" aria-label="Map controls">
          <button class="round-action" data-action="map-zoom-in" aria-label="Zoom in">+</button>
          <button class="round-action" data-action="map-zoom-out" aria-label="Zoom out">−</button>
          <button class="round-action" data-action="map-reset" aria-label="Reset map">${icon("focus")}</button>
          <span class="map-control-divider"></span>
          <button class="round-action ${state.darkMap ? "is-active" : ""}" data-action="toggle-map-theme" aria-label="Toggle map contrast">${icon("sun")}</button>
        </div>
        <div class="map-gesture-hint">${icon("hand")} Drag to pan · pinch or controls to zoom · tap a feature to inspect</div>
      </section>
      <aside class="map-inspector">
        ${renderDocumentMetrics(document, sites)}
        <section class="surface annotation-surface">
          <div class="surface-heading compact-heading"><div><span class="section-label">ANNOTATIONS</span><h2>Visible tracks</h2></div><button class="text-button" data-view="features">All ${document.features.length}</button></div>
          <div class="track-list">${visibleFeatures.slice(0, 6).map((item) => renderFeatureRow(item, document)).join("")}</div>
        </section>
      </aside>
    </div>
    <section class="surface restriction-strip">
      <div class="restriction-copy"><span class="tool-icon restriction-icon">${icon("scissors")}</span><div><span class="section-label">RESTRICTION OVERLAY</span><h2>${sites.length} visible cut sites across ${new Set(sites.map((site) => site.name)).size} enzymes</h2><p>${state.showSites ? "Tap a marker on the map or open Digest to inspect fragments." : "Markers are hidden; enable the overlay to inspect cut sites."}</p></div></div>
      <div class="restriction-actions"><button class="soft-action" data-action="toggle-sites">${state.showSites ? "Hide sites" : "Show sites"}</button><button class="primary-action" data-action="open-tool" data-tool="restriction">Analyze digest ${icon("arrow")}</button></div>
    </section>
  `;
}

function renderDocumentMetrics(document, sites) {
  return `<section class="surface stats-surface">
    <div class="surface-heading compact-heading"><div><span class="section-label">CONSTRUCT</span><h2>At a glance</h2></div><span class="topology-pill">${document.topology}</span></div>
    <div class="stat-grid">
      <div class="stat-card"><strong>${new Intl.NumberFormat().format(document.sequence.length)}</strong><span>base pairs</span></div>
      <div class="stat-card"><strong>${gcPercent(document.sequence)}%</strong><span>GC content</span></div>
      <div class="stat-card"><strong>${document.features.length}</strong><span>annotations</span></div>
      <div class="stat-card"><strong>${sites.length}</strong><span>cut sites</span></div>
    </div>
    <div class="origin-note"><span>${icon("pin")}</span><span>Origin at base 1 · revisions are local</span></div>
  </section>`;
}

function polarPoint(angle, radius, center = 250) {
  return { x: center + Math.cos(angle) * radius, y: center + Math.sin(angle) * radius };
}

function arcPath(start, end, radius, length, clockwise = true, center = 250) {
  const initial = -Math.PI / 2 + (start / length) * Math.PI * 2;
  let terminal = -Math.PI / 2 + (end / length) * Math.PI * 2;
  if (end <= start) terminal += Math.PI * 2;
  const from = polarPoint(initial, radius, center);
  const to = polarPoint(terminal, radius, center);
  const largeArc = terminal - initial > Math.PI ? 1 : 0;
  if (!clockwise) {
    return `M ${to.x.toFixed(2)} ${to.y.toFixed(2)} A ${radius} ${radius} 0 ${largeArc} 0 ${from.x.toFixed(2)} ${from.y.toFixed(2)}`;
  }
  return `M ${from.x.toFixed(2)} ${from.y.toFixed(2)} A ${radius} ${radius} 0 ${largeArc} 1 ${to.x.toFixed(2)} ${to.y.toFixed(2)}`;
}

function circularFeaturePath(item, document, index) {
  const radius = 160 + (index % 4) * 19;
  return featureSegments(item, document.sequence.length).map(([start, end]) => `<path d="${arcPath(start, end, radius, document.sequence.length, item.strand === 1)}" class="feature-arc ${state.activeFeatureId === item.id ? "is-selected" : ""}" style="stroke:${categoryColor(item)}" data-feature-id="${item.id}"/>`).join("");
}

function renderCircularMap(document, features, sites) {
  const tickAngles = Array.from({ length: 12 }, (_, index) => index * Math.PI * 2 / 12 - Math.PI / 2);
  return `<svg class="plasmid-map" id="plasmid-map" viewBox="0 0 500 500" role="img" aria-label="Circular plasmid map for ${escapeHtml(document.name)}">
      <defs>
        <filter id="map-shadow" x="-30%" y="-30%" width="160%" height="160%"><feDropShadow dx="0" dy="5" stdDeviation="7" flood-opacity=".13"/></filter>
        <linearGradient id="core-gradient" x1="0" y1="0" x2="1" y2="1"><stop stop-color="#f8fbf8"/><stop offset="1" stop-color="#e7f1eb"/></linearGradient>
      </defs>
      <circle class="map-halo" cx="250" cy="250" r="212"/>
      ${tickAngles.map((angle, index) => { const outer = polarPoint(angle, 207); const inner = polarPoint(angle, index % 3 === 0 ? 193 : 198); return `<line class="map-tick" x1="${outer.x}" y1="${outer.y}" x2="${inner.x}" y2="${inner.y}"/>`; }).join("")}
      <circle class="map-backbone" cx="250" cy="250" r="152"/>
      ${features.map((item, index) => circularFeaturePath(item, document, index)).join("")}
      ${state.showSites ? sites.map((site) => { const angle = -Math.PI / 2 + site.position / document.sequence.length * Math.PI * 2; const point = polarPoint(angle, 127); const outer = polarPoint(angle, 143); return `<g class="site-marker" data-site-position="${site.position}" data-site-name="${site.name}"><line x1="${point.x}" y1="${point.y}" x2="${outer.x}" y2="${outer.y}" style="stroke:${site.accent}"/><circle cx="${outer.x}" cy="${outer.y}" r="3.3" style="fill:${site.accent}"/></g>`; }).join("") : ""}
      <circle class="map-core" cx="250" cy="250" r="108" filter="url(#map-shadow)"/>
      <text x="250" y="235" class="map-name" text-anchor="middle">${escapeHtml(document.name)}</text>
      <text x="250" y="260" class="map-meta" text-anchor="middle">${new Intl.NumberFormat().format(document.sequence.length)} bp · ${document.topology}</text>
      <text x="250" y="283" class="map-revision" text-anchor="middle">OpenGene local revision r${document.revision}</text>
      <path class="origin-pointer" d="M 250 70 l -7 13 h 14 z"/>
    </svg>`;
}

function renderLinearMap(document, features, sites) {
  const featureBlocks = features.map((item, index) => featureSegments(item, document.sequence.length).map(([start, end]) => {
    const left = start / document.sequence.length * 100;
    const width = Math.max(1.2, (end - start) / document.sequence.length * 100);
    const row = index % 4;
    return `<button class="linear-feature ${state.activeFeatureId === item.id ? "is-selected" : ""}" data-feature-id="${item.id}" style="--left:${left}%;--width:${width}%;--row:${row};--feature-color:${categoryColor(item)}"><span>${escapeHtml(item.name)}</span></button>`;
  }).join("")).join("");
  return `<div class="linear-map" id="linear-map">
    <div class="linear-ruler"><span>1</span><span>${formatPosition(Math.round(document.sequence.length * .25))}</span><span>${formatPosition(Math.round(document.sequence.length * .5))}</span><span>${formatPosition(Math.round(document.sequence.length * .75))}</span><span>${new Intl.NumberFormat().format(document.sequence.length)}</span></div>
    <div class="linear-backbone"></div>
    ${featureBlocks}
    ${state.showSites ? sites.map((site) => `<button class="linear-site" data-site-position="${site.position}" data-site-name="${site.name}" style="--site:${site.position / document.sequence.length * 100}%;--site-color:${site.accent}" title="${site.name} at ${formatPosition(site.position)}"></button>`).join("") : ""}
    <div class="linear-center-note">${escapeHtml(document.name)} · ${document.topology} DNA</div>
  </div>`;
}

function renderFeatureRow(item, document) {
  return `<button class="feature-row ${state.activeFeatureId === item.id ? "is-active" : ""}" data-feature-id="${item.id}">
      <span class="feature-swatch" style="--feature-color:${categoryColor(item)}"></span>
      <span class="feature-row-main"><strong>${escapeHtml(item.name)}</strong><small>${escapeHtml(item.type.replaceAll("_", " "))} · ${formatBp(featureLength(item, document))}</small></span>
      <span class="feature-arrow">${item.strand === 1 ? "→" : "←"}</span>
    </button>`;
}

function renderSequenceView(document) {
  const selectionStart = state.selection.start;
  const activeFeature = document.features.find((item) => item.id === state.activeFeatureId);
  return `<div class="sequence-layout">
      <section class="surface sequence-surface">
        <div class="surface-heading sequence-heading"><div><span class="section-label">BASE-LEVEL EDITOR</span><h2>Sequence canvas</h2></div><div class="sequence-header-actions"><span class="mono-chip">${new Intl.NumberFormat().format(document.sequence.length)} nt</span><button class="soft-action" data-action="edit-sequence">${icon("edit")} Edit selection</button></div></div>
        <div class="sequence-guide"><span>Tap a base to move the cursor.</span><span>Drag across bases to select.</span><span>10 bases per group · 60 per row</span></div>
        <div class="sequence-scroll" id="sequence-scroll" tabindex="0" aria-label="DNA sequence editor">${renderSequenceRows(document)}</div>
      </section>
      <aside class="sequence-inspector">
        <section class="surface selection-surface">
          <div class="surface-heading compact-heading"><div><span class="section-label">SELECTION</span><h2>${selectionLength(document) ? formatBp(selectionLength(document)) : "Cursor"}</h2></div><button class="icon-action" data-action="copy-selection" aria-label="Copy selection">${icon("copy")}</button></div>
          <div class="selection-coordinates"><span>FROM <b>${formatPosition(selectionStart)}</b></span><span>TO <b>${formatPosition(state.selection.end || selectionStart)}</b></span></div>
          <code class="sequence-preview">${escapeHtml(selectionSequence(document).slice(0, 72) || document.sequence.slice(Math.max(0, selectionStart - 10), selectionStart + 30))}</code>
          <div class="selection-buttons"><button class="soft-action" data-action="reverse-complement-selection">${icon("reverse")} Reverse complement</button><button class="soft-action" data-action="translate-selection">${icon("translate")} Translate</button></div>
        </section>
        <section class="surface feature-inspector">
          <div class="surface-heading compact-heading"><div><span class="section-label">AT CURSOR</span><h2>${activeFeature ? escapeHtml(activeFeature.name) : "No annotation"}</h2></div>${activeFeature ? `<button class="icon-action danger-action" data-action="delete-feature" data-feature-id="${activeFeature.id}" aria-label="Delete annotation">${icon("trash")}</button>` : ""}</div>
          ${activeFeature ? `<div class="feature-detail-color" style="--feature-color:${categoryColor(activeFeature)}"></div><dl class="feature-dl"><div><dt>Type</dt><dd>${escapeHtml(activeFeature.type.replaceAll("_", " "))}</dd></div><div><dt>Location</dt><dd>${locationLabel(activeFeature, document)}</dd></div><div><dt>Length</dt><dd>${formatBp(featureLength(activeFeature, document))}</dd></div></dl><p>${escapeHtml(activeFeature.qualifiers.note || activeFeature.qualifiers.gene || "No additional note.")}</p>` : `<p class="muted-copy">Select an annotation from the map or feature list to inspect its location and qualifiers.</p>`}
        </section>
      </aside>
    </div>`;
}

function renderSequenceRows(document) {
  const lines = [];
  const size = 60;
  for (let start = 0; start < document.sequence.length; start += size) {
    const row = document.sequence.slice(start, start + size);
    const chunks = [];
    for (let groupStart = 0; groupStart < row.length; groupStart += 10) {
      const group = row.slice(groupStart, groupStart + 10);
      chunks.push(`<span class="base-group">${[...group].map((base, index) => {
        const position = start + groupStart + index;
        const isSelected = isPositionSelected(position, document);
        const featureAtBase = document.features.find((item) => featureContains(item, position, document));
        return `<button class="base ${isSelected ? "is-selected" : ""} ${featureAtBase ? "has-feature" : ""}" data-base-position="${position}" aria-label="Base ${formatPosition(position)} ${base}">${base}</button>`;
      }).join("")}</span>`);
    }
    lines.push(`<div class="sequence-row" data-sequence-row="${start}"><span class="base-index">${String(start + 1).padStart(6)}</span><div class="base-groups">${chunks.join("")}</div><span class="base-index end-index">${String(Math.min(start + size, document.sequence.length)).padStart(6)}</span></div>`);
  }
  return lines.join("");
}

function isPositionSelected(position, document) {
  const { start, end } = state.selection;
  if (start === end) return position === start;
  if (end > start) return position >= start && position < end;
  return position >= start || position < end;
}

function featureContains(item, position, document) {
  return featureSegments(item, document.sequence.length).some(([start, end]) => position >= start && position < end);
}

function renderFeaturesView(document) {
  const selected = document.features.find((item) => item.id === state.activeFeatureId);
  return `<div class="features-layout">
      <section class="surface feature-table-surface">
        <div class="surface-heading features-heading"><div><span class="section-label">ANNOTATION LIBRARY</span><h2>${document.features.length} features on this construct</h2></div><button class="primary-action" data-action="add-feature">${icon("plus")} Add feature</button></div>
        <label class="search-box">${icon("search")}<input id="feature-search" value="${escapeHtml(state.search)}" placeholder="Filter by name or type" autocomplete="off"/><kbd>⌘ K</kbd></label>
        <div class="feature-table" role="table"><div class="feature-table-head" role="row"><span>Feature</span><span>Type</span><span>Coordinates</span><span>Length</span><span></span></div><div id="feature-table-body">${renderFeatureTableRows(document, selected)}</div></div>
      </section>
      <aside class="feature-detail-panel">${renderFeatureDetail(selected, document)}</aside>
    </div>`;
}

function renderFeatureTableRows(document, selected) {
  const filtered = document.features.filter((item) => `${item.name} ${item.type}`.toLowerCase().includes(state.search.toLowerCase()));
  return filtered.map((item) => `<div class="feature-table-row ${selected?.id === item.id ? "is-selected" : ""}" role="row"><button data-feature-id="${item.id}" class="feature-name-cell"><span class="feature-swatch" style="--feature-color:${categoryColor(item)}"></span><span><strong>${escapeHtml(item.name)}</strong><small>${item.strand === 1 ? "Forward strand" : "Reverse strand"}</small></span></button><span class="type-cell">${escapeHtml(item.type.replaceAll("_", " "))}</span><span class="coords-cell">${locationLabel(item, document)}</span><span>${formatBp(featureLength(item, document))}</span><button class="row-more" data-action="feature-menu" data-feature-id="${item.id}" aria-label="Feature actions">${icon("dots")}</button></div>`).join("") || `<div class="no-results">No annotations match “${escapeHtml(state.search)}”.</div>`;
}

function renderFeatureDetail(item, document) {
  if (!item) return `<section class="surface empty-detail"><div class="empty-orb small-orb">${icon("tag")}</div><h2>Inspect an annotation</h2><p>Tap a row to see qualifiers, coordinates, and editing actions.</p></section>`;
  return `<section class="surface detail-surface"><div class="feature-detail-top"><span class="feature-large-swatch" style="--feature-color:${categoryColor(item)}"></span><div><span class="section-label">${escapeHtml(item.type.replaceAll("_", " "))}</span><h2>${escapeHtml(item.name)}</h2></div><button class="icon-action" data-action="delete-feature" data-feature-id="${item.id}" aria-label="Remove feature">${icon("trash")}</button></div><div class="detail-location"><span>${icon("pin")}</span><div><strong>${locationLabel(item, document)}</strong><small>${formatBp(featureLength(item, document))} · ${item.strand === 1 ? "forward" : "reverse"} strand</small></div></div><dl class="feature-dl">${Object.entries(item.qualifiers || {}).filter(([key]) => key !== "wrapsOrigin").map(([key, value]) => `<div><dt>${escapeHtml(key)}</dt><dd>${escapeHtml(value)}</dd></div>`).join("") || `<div><dt>Note</dt><dd>No qualifiers supplied</dd></div>`}</dl><div class="detail-actions"><button class="soft-action" data-action="focus-feature" data-feature-id="${item.id}">${icon("focus")} Focus on map</button><button class="soft-action" data-action="copy-feature-sequence" data-feature-id="${item.id}">${icon("copy")} Copy bases</button></div></section>`;
}

function translateBrowserSequence(sequence) {
  let protein = "";
  for (let index = 0; index + 2 < sequence.length; index += 3) protein += GENETIC_CODE[sequence.slice(index, index + 3)] || "X";
  return protein;
}

function browserOrfCandidates(document, minimumProteinLength = 30) {
  const sequence = document.sequence;
  const length = sequence.length;
  if (length < 6) return [];
  const circular = document.topology === "circular";
  const candidates = [];
  const scan = (oriented, strand) => {
    const searchable = circular ? `${oriented}${oriented}` : oriented;
    for (let frame = 0; frame < 3; frame += 1) {
      for (let start = frame; start + 2 < length; start += 3) {
        if (searchable.slice(start, start + 3) !== "ATG") continue;
        const limit = circular ? start + length : length;
        for (let stopStart = start + 3; stopStart + 2 < limit; stopStart += 3) {
          const stop = searchable.slice(stopStart, stopStart + 3);
          if (!STOP_CODONS.has(stop)) continue;
          const end = stopStart + 3;
          const codingSequence = searchable.slice(start, stopStart);
          const protein = translateBrowserSequence(codingSequence);
          if (protein.length < minimumProteinLength) break;
          let featureStart;
          let featureEnd;
          let wrapsOrigin = false;
          if (strand === 1) {
            featureStart = start;
            featureEnd = end <= length ? end : end % length;
            wrapsOrigin = end > length;
          } else {
            const endModulo = end % length;
            featureStart = endModulo === 0 ? 0 : length - endModulo;
            featureEnd = start === 0 ? length : length - start;
            wrapsOrigin = end > length;
          }
          candidates.push({
            key: `${strand}:${start}:${end}`,
            start: featureStart,
            end: featureEnd,
            strand,
            wrapsOrigin,
            frame: strand === 1 ? frame + 1 : -(frame + 1),
            protein,
            nucleotideLength: end - start,
            stopCodon: stop,
          });
          break;
        }
      }
    }
  };
  scan(sequence, 1);
  scan(reverseComplement(sequence), -1);
  return candidates.sort((left, right) => right.protein.length - left.protein.length || left.start - right.start || right.strand - left.strand);
}

function addBrowserOrf(index) {
  const document = activeDocument();
  const candidate = browserOrfCandidates(document)[index];
  if (!candidate) return showToast("That ORF suggestion is no longer available");
  if (document.features.some((item) => item.qualifiers?.sourceOrfKey === candidate.key)) return showToast("This ORF is already annotated");
  const item = feature("cds", `ORF ${candidate.frame > 0 ? "+" : ""}${candidate.frame} · ${candidate.protein.length} aa`, candidate.start, candidate.end, candidate.strand, {
    sourceOrfKey: candidate.key,
    translation: candidate.protein,
    start_codon: "ATG",
    terminal_stop_codon: candidate.stopCodon,
    wrapsOrigin: candidate.wrapsOrigin,
    note: "Local six-frame ORF suggestion — review before treating it as an annotation.",
  });
  commit(`Added ${item.name}`, () => {
    document.features.push(item);
    state.activeFeatureId = item.id;
    state.selection = { start: item.start, end: item.end };
  }, `ORF scan suggestion · ${formatBp(candidate.nucleotideLength)} · frame ${candidate.frame > 0 ? "+" : ""}${candidate.frame}`);
}

function annotateBrowserOrfs() {
  const document = activeDocument();
  const candidates = browserOrfCandidates(document).filter((candidate) => !document.features.some((item) => item.qualifiers?.sourceOrfKey === candidate.key));
  if (!candidates.length) return showToast("No new 30 aa ORF suggestions were found");
  const selected = candidates.slice(0, 24);
  commit(`Added ${selected.length} ORF suggestion${selected.length === 1 ? "" : "s"}`, () => {
    for (const candidate of selected) {
      document.features.push(feature("cds", `ORF ${candidate.frame > 0 ? "+" : ""}${candidate.frame} · ${candidate.protein.length} aa`, candidate.start, candidate.end, candidate.strand, {
        sourceOrfKey: candidate.key,
        translation: candidate.protein,
        start_codon: "ATG",
        terminal_stop_codon: candidate.stopCodon,
        wrapsOrigin: candidate.wrapsOrigin,
        note: "Local six-frame ORF suggestion — review before treating it as an annotation.",
      }));
    }
  }, `Six-frame ORF scan · ${selected.length} evidence-bearing CDS suggestions`);
}

function autoAnnotate() {
  const document = activeDocument();
  if (!document) return;
  const circular = document.topology === "circular";
  const searchable = circular ? document.sequence + document.sequence.slice(0, document.sequence.length - 1) : document.sequence;
  let added = 0;

  commit(`Auto-annotate features`, () => {
    COMMON_FEATURES_DB.forEach((dbFeature) => {
      const dbSeq = normalizeSequence(dbFeature.sequence);
      if (dbSeq.length < 10) return;
      const orientations = [
        { seq: dbSeq, strand: 1 },
        { seq: reverseComplement(dbSeq), strand: -1 }
      ];

      orientations.forEach(orientation => {
        let index = searchable.indexOf(orientation.seq);
        while (index !== -1) {
          if (index < document.sequence.length) {
            let start = index;
            let end = (index + orientation.seq.length) % document.sequence.length;
            if (end === 0) end = document.sequence.length;
            const exists = document.features.some(f => f.name === dbFeature.name && f.start === start && f.end === end && f.strand === orientation.strand);

            if (!exists) {
              const wrapsOrigin = (index + orientation.seq.length) > document.sequence.length;
              const f = feature(dbFeature.type, dbFeature.name, start, end, orientation.strand, { note: "Auto-annotated common feature", wrapsOrigin });
              document.features.push(f);
              added++;
            }
          }
          index = searchable.indexOf(orientation.seq, index + 1);
        }
      });
    });
  });

  if (added) {
    showToast(`Added ${added} annotations`);
    render();
  } else {
    showToast("No new common features found");
  }
}

function renderAnnotationTool(document) {
  const candidates = browserOrfCandidates(document);
  const newCount = candidates.filter((candidate) => !document.features.some((item) => item.qualifiers?.sourceOrfKey === candidate.key)).length;
  const addableCount = Math.min(newCount, 24);
  return `<div class="tool-heading"><div><span class="section-label">ANNOTATION &amp; CODING</span><h2>Inspect coding potential</h2><p>Local six-frame scan shows complete ATG-to-stop ORF suggestions of at least 30 amino acids. Suggestions are never silently accepted.</p></div><div class="tool-heading-actions"><button class="soft-action" data-action="auto-annotate">${icon("spark")} Auto-annotate</button><button class="soft-action" data-action="open-codon-optimize">${icon("translate")} Reverse translate</button><button class="primary-action" data-action="annotate-all-orfs" ${newCount ? "" : "disabled"}>${icon("plus")} Add ${addableCount || ""} suggestion${addableCount === 1 ? "" : "s"}</button></div></div>
    <div class="annotation-summary"><section><span class="section-label">ORF SUGGESTIONS</span><strong>${candidates.length}</strong><small>${newCount} not yet annotated</small></section><section><span class="section-label">CODING FEATURES</span><strong>${document.features.filter((item) => ["cds", "gene"].includes(item.type)).length}</strong><small>existing annotations</small></section><section><span class="section-label">TOPOLOGY</span><strong>${document.topology}</strong><small>${document.topology === "circular" ? "origin-spanning candidates included" : "linear six-frame scan"}</small></section></div>
    ${candidates.length ? `<div class="orf-list">${candidates.slice(0, 10).map((candidate, index) => `<article class="orf-row"><span class="orf-frame ${candidate.strand === 1 ? "forward" : "reverse"}">${candidate.frame > 0 ? "+" : ""}${candidate.frame}</span><div><strong>ORF · ${candidate.protein.length} aa</strong><small>${formatPosition(candidate.start)}${candidate.wrapsOrigin ? " → " : " – "}${formatPosition(candidate.wrapsOrigin ? candidate.end : Math.max(candidate.start, candidate.end - 1))}${candidate.wrapsOrigin ? " (wrap)" : ""} · ${formatBp(candidate.nucleotideLength)} · stop ${candidate.stopCodon}</small><code>${escapeHtml(candidate.protein.slice(0, 54))}${candidate.protein.length > 54 ? "…" : ""}</code></div><button class="soft-action" data-action="add-orf" data-orf-index="${index}">${icon("plus")} Add</button></article>`).join("")}</div>${candidates.length > 10 ? `<p class="tool-footnote">Showing the ten longest candidates. The portable Python API exposes all ${candidates.length} suggestions with structured evidence.</p>` : ""}` : `<section class="no-read-state"><div class="empty-orb">${icon("translate")}</div><h3>No complete ORFs at this threshold</h3><p>Try a different sequence or use the portable API to scan shorter and partial coding regions.</p><button class="primary-action" data-action="open-codon-optimize">${icon("translate")} Reverse translate protein</button></section>`}
    <section class="annotation-note">${icon("info")} <span>Codon design uses a transparent generic preference table in this browser. Choose organism-specific tables and run independent scientific review before ordering DNA.</span></section>`;
}

function renderToolsView(document) {
  const tools = [
    ["restriction", "Restriction digest", "Scan enzymes, inspect fragments, and preview a gel", "scissors"],
    ["primers", "Primer & PCR", "Design primers and simulate an amplicon", "primer"],
    ["annotations", "Annotation & coding", "Find ORFs, translate, and make coding candidates", "translate"],
    ["cloning", "Assembly planning", "Create a reproducible, non-mutating operation record", "assembly"],
    ["verification", "Verify reads", "Align reads, review variants, and build a consensus", "verify"],
    ["gel", "Virtual gel", "Compare digestion and PCR bands", "gel"],
  ];
  return `<div class="tools-layout"><aside class="tool-rail"><span class="section-label">ANALYSIS TOOLS</span>${tools.map(([id, title, description, glyph]) => `<button data-action="open-tool" data-tool="${id}" class="tool-nav ${state.activeTool === id ? "is-active" : ""}"><span class="tool-nav-icon">${icon(glyph)}</span><span><strong>${title}</strong><small>${description}</small></span></button>`).join("")}</aside><section class="tool-workspace surface">${renderToolPanel(document)}</section></div>`;
}

function renderToolPanel(document) {
  if (state.activeTool === "restriction") return renderRestrictionTool(document);
  if (state.activeTool === "primers") return renderPrimerTool(document);
  if (state.activeTool === "annotations") return renderAnnotationTool(document);
  if (state.activeTool === "cloning") return renderCloningTool(document);
  if (state.activeTool === "verification") return renderVerificationTool(document);
  return renderGelTool(document);
}

function renderRestrictionTool(document) {
  const sites = allSites(document);
  const selectedEnzymes = ENZYMES.slice(0, 3);
  const selectedSites = sites.filter((site) => selectedEnzymes.some((enzyme) => enzyme.name === site.name));
  const fragments = digestFragments(document, selectedSites).sort((a, b) => b.length - a.length);
  return `<div class="tool-heading"><div><span class="section-label">RESTRICTION ANALYSIS</span><h2>Digest planner</h2><p>Independent recognition-site scan. Select enzymes to preview cut positions and fragment sizes.</p></div><button class="primary-action" data-action="run-digest">${icon("play")} Run digest</button></div><div class="enzyme-picker">${ENZYMES.map((enzyme) => `<button class="enzyme-chip ${selectedEnzymes.some((item) => item.name === enzyme.name) ? "is-selected" : ""}" data-action="enzyme-toggle" data-enzyme="${enzyme.name}"><span style="--enzyme-color:${enzyme.accent}"></span><strong>${enzyme.name}</strong><small>${enzyme.site}</small></button>`).join("")}</div><div class="digest-layout"><section class="digest-results"><div class="digest-summary"><div><span class="section-label">CUT SUMMARY</span><h3>${selectedSites.length} cuts · ${plural(fragments.length, "fragment")}</h3></div><span class="topology-pill">${document.topology}</span></div><div class="fragment-list">${fragments.map((fragment, index) => `<button class="fragment-row" data-action="select-fragment" data-start="${fragment.start}" data-end="${fragment.end}"><span class="fragment-number">${index + 1}</span><span class="fragment-scale"><i style="--fragment-width:${Math.max(5, fragment.length / document.sequence.length * 100)}%"></i></span><strong>${formatBp(fragment.length)}</strong><small>${formatPosition(fragment.start)} → ${formatPosition(fragment.end)}</small></button>`).join("")}</div></section><section class="digest-map-card">${renderDigestMiniMap(document, selectedSites)}<div class="mini-map-legend">${selectedEnzymes.map((enzyme) => `<span><i style="--enzyme-color:${enzyme.accent}"></i>${enzyme.name}</span>`).join("")}</div></section></div>`;
}

function renderDigestMiniMap(document, sites) {
  return `<svg class="digest-mini-map" viewBox="0 0 260 260" aria-label="Restriction digest site overview"><circle cx="130" cy="130" r="78" class="digest-ring"/>${sites.map((site) => { const angle = -Math.PI / 2 + site.position / document.sequence.length * Math.PI * 2; const one = polarPoint(angle, 65, 130); const two = polarPoint(angle, 91, 130); return `<g><line x1="${one.x}" y1="${one.y}" x2="${two.x}" y2="${two.y}" style="stroke:${site.accent}"/><circle cx="${two.x}" cy="${two.y}" r="3.5" style="fill:${site.accent}"/></g>`; }).join("")}<text x="130" y="126" text-anchor="middle">${sites.length} CUTS</text><text x="130" y="146" text-anchor="middle">${document.name}</text></svg>`;
}

function renderPrimerTool(document) {
  const primers = primerCandidates(document);
  return `<div class="tool-heading"><div><span class="section-label">PRIMER & PCR</span><h2>Design around your selection</h2><p>Suggestions use local, explainable length, GC, and melting-temperature checks. Validate experimentally.</p></div><button class="primary-action" data-action="design-primers">${icon("spark")} Refresh suggestions</button></div><div class="primer-context"><span>${icon("cursor")}</span><div><strong>${selectionLength(document) ? `${formatBp(selectionLength(document))} selection` : "No range selected"}</strong><small>Choose bases in Sequence view to anchor an amplicon.</small></div><button class="text-button" data-view="sequence">Choose range</button></div><div class="primer-grid">${primers.map((primer) => `<section class="primer-card"><div class="primer-card-top"><span class="primer-direction ${primer.direction.toLowerCase()}">${primer.direction}</span><button class="icon-action" data-action="copy-primer" data-primer="${primer.sequence}" aria-label="Copy ${primer.name}">${icon("copy")}</button></div><h3>${primer.name}</h3><code>${primer.sequence}</code><div class="primer-metrics"><span><b>${primer.tm}°C</b> Tm</span><span><b>${primer.gc}%</b> GC</span><span><b>${primer.sequence.length}</b> nt</span></div></section>`).join("")}</div><section class="pcr-result"><div><span class="section-label">VIRTUAL PCR</span><h3>Predicted product: ${formatBp(primers[0]?.product || 0)}</h3><p>Forward and reverse candidates are within a 3°C Tm window.</p></div><div class="pcr-status"><span class="pass-dot"></span> Candidate pair passes basic screen</div><button class="soft-action" data-action="annotate-primers">${icon("tag")} Add to map</button></section>`;
}

function renderCloningTool(document) {
  const methods = [
    ["Gibson assembly", "Overlap-driven, scarless assembly", "gibson"],
    ["Golden Gate", "Type IIS modular assembly", "golden"],
    ["Restriction + ligation", "Classic digest and ligation", "ligation"],
    ["Gateway transfer", "att-site recombination workflow", "gateway"],
    ["In-Fusion / overlap", "Homology-guided insert planning", "infusion"],
    ["TA / blunt cloning", "End-compatible insertion planning", "ta"],
    ["Oligo annealing", "Create a duplex insert from paired oligos", "oligo"],
    ["Synthetic insert", "Add a defined fragment at the selection", "synthetic"],
  ];
  return `<div class="tool-heading"><div><span class="section-label">ASSEMBLY PLANNING</span><h2>Build a reproducible operation</h2><p>Each workflow records inputs, parameters, warnings, and reproducibility metadata before creating a candidate product.</p></div><button class="primary-action" data-action="create-assembly">${icon("plus")} New operation</button></div><div class="assembly-steps"><div class="assembly-step is-complete"><span>1</span><div><strong>Choose inputs</strong><small>${escapeHtml(document.name)} is the active backbone</small></div><i>${icon("check")}</i></div><div class="assembly-step is-current"><span>2</span><div><strong>Select an assembly method</strong><small>Plans are deterministic and never overwrite the source construct.</small></div><i>${icon("arrow")}</i></div><div class="assembly-step"><span>3</span><div><strong>Review candidate product</strong><small>Inspect sequence, warnings, and the saved provenance node.</small></div></div></div><div class="assembly-methods expanded-methods">${methods.map(([title, detail, method]) => `<button class="assembly-method" data-action="create-assembly" data-method="${method}"><span class="method-icon">${icon("assembly")}</span><span><strong>${title}</strong><small>${detail}</small></span>${icon("arrow")}</button>`).join("")}</div><div class="provenance-note">${icon("history")} The native and automation engines use the same immutable input/output contract. Review wet-lab parameters before use.</div>`;
}

function renderGelTool(document) {
  const fragments = digestFragments(document, allSites(document).filter((site) => ["EcoRI", "BamHI", "HindIII"].includes(site.name))).sort((a, b) => b.length - a.length);
  return `<div class="tool-heading"><div><span class="section-label">VIRTUAL GEL</span><h2>Fragment visualization</h2><p>A deterministic size-to-position preview; it is not a laboratory-grade migration prediction.</p></div><button class="primary-action" data-action="run-gel">${icon("play")} Render gel</button></div><div class="gel-settings"><label>Gel percentage<select id="gel-percent"><option>0.8%</option><option selected>1.0%</option><option>1.5%</option><option>2.0%</option></select></label><label>Lane width<select id="gel-width"><option>standard</option><option>wide</option></select></label><span class="gel-help">${icon("info")} Fragment intensity tracks relative mass.</span></div><div class="gel-stage"><div class="gel-labels"><span>10 kb</span><span>3 kb</span><span>1 kb</span><span>500 bp</span><span>100 bp</span></div><div class="gel-image"><div class="gel-lane ladder-lane"><span class="lane-title">Ladder</span>${[10000, 6000, 3000, 2000, 1000, 700, 500, 300, 100].map((length) => gelBand(length, 10000, "ladder")).join("")}</div><div class="gel-lane"><span class="lane-title">Digest</span>${fragments.map((fragment) => gelBand(fragment.length, document.sequence.length, "digest")).join("")}</div><div class="gel-lane"><span class="lane-title">PCR</span>${gelBand(Math.max(120, selectionLength(document) || 650), document.sequence.length, "pcr")}</div></div><div class="gel-caption">EcoRI + BamHI + HindIII digest · ${plural(fragments.length, "fragment")}</div></div>`;
}

function renderVerificationTool(document) {
  const summary = verificationSummary(document);
  const bins = coverageBins(summary);
  const status = confidenceLabel(summary.meanIdentity);
  const mapped = summary.alignments.filter((entry) => entry.alignment);
  return `<div class="tool-heading"><div><span class="section-label">SEQUENCE VERIFICATION</span><h2>Align reads to ${escapeHtml(document.name)}</h2><p>Local, explainable ungapped read mapping for quick review. Use the portable Python alignment workflow for indel-aware scientific verification.</p></div><div class="tool-heading-actions"><button class="soft-action" data-action="add-read">${icon("upload")} Add read</button><button class="primary-action" data-action="run-verification">${icon("play")} Align ${state.reads.length} reads</button></div></div>
    ${state.reads.length ? `<div class="verification-overview"><section class="verification-metric"><span class="section-label">MEAN IDENTITY</span><strong>${(summary.meanIdentity * 100).toFixed(1)}%</strong><small>${status}</small></section><section class="verification-metric"><span class="section-label">REFERENCE COVERED</span><strong>${document.sequence.length ? (summary.coveredBases / document.sequence.length * 100).toFixed(1) : 0}%</strong><small>${new Intl.NumberFormat().format(summary.coveredBases)} / ${new Intl.NumberFormat().format(document.sequence.length)} bp</small></section><section class="verification-metric"><span class="section-label">CALLED VARIANTS</span><strong>${summary.variants.length}</strong><small>${summary.variants.length ? "review before acceptance" : "no confident substitutions"}</small></section></div>
    <section class="coverage-card"><div class="coverage-heading"><div><span class="section-label">COVERAGE</span><h3>Read depth across reference</h3></div><span class="quality-pill">${mapped.length}/${state.reads.length} reads mapped</span></div><div class="coverage-plot" aria-label="Reference coverage plot">${bins.map((depth, index) => `<i class="coverage-bar ${depth ? "has-depth" : ""}" style="--coverage-height:${Math.min(100, 16 + depth * 34)}%" title="Region ${index + 1}: ${depth.toFixed(1)}×"></i>`).join("")}</div><div class="coverage-axis"><span>1</span><span>${new Intl.NumberFormat().format(Math.round(document.sequence.length / 2))}</span><span>${new Intl.NumberFormat().format(document.sequence.length)}</span></div></section>
    <div class="verification-grid"><section class="read-table-card"><div class="verification-card-heading"><div><span class="section-label">READS</span><h3>${state.reads.length} local sequence read${state.reads.length === 1 ? "" : "s"}</h3></div><button class="text-button" data-action="add-read">Add another</button></div><div class="read-table">${summary.alignments.map(({ read, alignment }) => `<button class="read-row" data-action="inspect-read" data-read-id="${read.id}"><span class="read-status ${alignment?.identity >= .97 ? "pass" : "review"}">${alignment ? icon(alignment.identity >= .97 ? "check" : "info") : icon("close")}</span><span class="read-name"><strong>${escapeHtml(read.name)}</strong><small>${new Intl.NumberFormat().format(read.sequence.length)} nt · ${escapeHtml(read.source || "manual")}</small></span><span class="read-map">${alignment ? `${formatPosition(alignment.start)} · ${alignment.strand}<b>${(alignment.identity * 100).toFixed(1)}%</b>` : "unmapped"}</span>${icon("arrow")}</button>`).join("")}</div></section><section class="variant-card"><div class="verification-card-heading"><div><span class="section-label">VARIANT REVIEW</span><h3>${summary.variants.length ? `${summary.variants.length} confident substitution${summary.variants.length === 1 ? "" : "s"}` : "No confident substitutions"}</h3></div><button class="text-button" data-action="create-consensus" ${mapped.length ? "" : "disabled"}>Create consensus</button></div>${summary.variants.length ? `<div class="variant-list">${summary.variants.slice(0, 6).map((variant) => `<button class="variant-row" data-action="focus-variant" data-position="${variant.position}"><span>${formatPosition(variant.position)}</span><strong>${variant.reference} → ${variant.alternate}</strong><small>${variant.support}/${variant.depth} reads</small></button>`).join("")}</div>` : `<div class="empty-verification">${icon("check")} <span>All confident calls match the active reference.</span></div>`}</section></div>
    <section class="verification-note">${icon("info")} <span><strong>Trace files:</strong> paste nucleotide reads or import FASTQ here. The Python sequencing package exposes an optional real-AB1 parser boundary; this browser UI does not fabricate chromatograms or indel calls.</span></section>` : `<section class="no-read-state"><div class="empty-orb">${icon("verify")}</div><h3>Add sequence reads for verification</h3><p>Paste a FASTA-like read sequence, then map it locally against the active construct. Imported traces remain explicit, never inferred.</p><button class="primary-action" data-action="add-read">${icon("upload")} Add first read</button></section>`}`;
}

function gelBand(length, max, type) {
  const position = Math.min(92, Math.max(5, 8 + (Math.log10(max) - Math.log10(Math.max(80, length))) / Math.max(1, Math.log10(max) - 2) * 78));
  const opacity = Math.min(.95, Math.max(.35, Math.sqrt(length / max) * 1.8));
  return `<i class="gel-band ${type}" style="--gel-y:${position}%;--gel-opacity:${opacity}"><b>${length >= 1000 ? `${(length / 1000).toFixed(length % 1000 ? 1 : 0)}k` : length}</b></i>`;
}

function renderHistoryView(document) {
  const entries = state.history;
  return `<div class="history-layout"><section class="surface history-canvas"><div class="surface-heading"><div><span class="section-label">CONSTRUCT PROVENANCE</span><h2>Reproducible history</h2><p>Every edit and operation is an explicit node in the local project DAG.</p></div><button class="soft-action" data-action="export-manifest">${icon("download")} Export manifest</button></div><div class="history-graph">${entries.slice().reverse().map((entry, index) => `<div class="history-node ${entry.type}" data-history-id="${entry.id}"><span class="history-connector"></span><div class="history-node-icon">${icon(entry.type === "edit" ? "edit" : entry.type === "cloning" ? "assembly" : entry.type === "import" ? "upload" : "tag")}</div><div class="history-node-body"><span class="history-time">${formatRelative(entry.timestamp)}</span><strong>${escapeHtml(entry.title)}</strong><small>${escapeHtml(entry.detail || "No operation detail")}</small>${entry.warnings?.length ? `<em>${entry.warnings.length} warning${entry.warnings.length > 1 ? "s" : ""}</em>` : ""}</div><button class="icon-action" data-action="history-detail" data-history-id="${entry.id}" aria-label="Inspect operation">${icon("arrow")}</button></div>`).join("")}</div></section><aside class="history-side"><section class="surface reproducibility-card"><span class="repro-icon">${icon("check")}</span><span class="section-label">REPRODUCIBILITY</span><h2>Local provenance intact</h2><p>${entries.length} operation records · source sequences are retained in the project snapshot.</p><div><span class="pass-dot"></span> Inputs and parameters captured</div></section><section class="surface hash-card"><span class="section-label">ACTIVE SEQUENCE HASH</span><code>${simpleHash(document.sequence)}</code><small>Deterministic local checksum · r${document.revision}</small></section></aside></div>`;
}

function formatRelative(timestamp) {
  const minutes = Math.max(0, Math.floor((Date.now() - new Date(timestamp).getTime()) / 60000));
  if (minutes < 1) return "just now";
  if (minutes < 60) return `${minutes} min ago`;
  if (minutes < 24 * 60) return `${Math.floor(minutes / 60)} hr ago`;
  return new Date(timestamp).toLocaleDateString(undefined, { month: "short", day: "numeric" });
}

function simpleHash(value) {
  let hash = 2166136261;
  for (let index = 0; index < value.length; index += 1) {
    hash ^= value.charCodeAt(index);
    hash = Math.imul(hash, 16777619);
  }
  return `fnv1a-${(hash >>> 0).toString(16).padStart(8, "0")}-${value.length}`;
}

function renderDocumentTree() {
  const groups = [
    ["Circular constructs", state.documents.filter((item) => item.topology === "circular")],
    ["Linear fragments", state.documents.filter((item) => item.topology !== "circular")],
  ].filter(([, documents]) => documents.length);
  return groups.map(([title, documents]) => `
    <div class="tree-row tree-row--group" aria-label="${escapeHtml(title)}">
      <span class="disclosure"><svg><use href="#i-chevron-down"></use></svg></span>
      <span class="tree-row__type tree-row__type--folder"><svg><use href="#i-folder"></use></svg></span>
      <span class="tree-row__label">${escapeHtml(title)}</span>
      <span class="tree-row__count">${documents.length}</span>
    </div>
    <div class="tree-children">
      ${documents.map((item) => `<button class="tree-row ${item.id === state.activeId ? "is-selected" : ""}" type="button" ${item.id === state.activeId ? "aria-current=\"page\"" : ""} data-action="open-document" data-document-id="${item.id}">
        <span class="tree-row__type tree-row__type--${item.topology === "circular" ? "plasmid" : "fragment"}"><svg><use href="#i-dna"></use></svg></span>
        <span class="tree-row__label"><b>${escapeHtml(item.name)}</b><small>${new Intl.NumberFormat().format(item.sequence.length)} bp · ${escapeHtml(item.topology)}${item.metadata?.consensus ? " · consensus" : item.metadata?.candidate ? " · candidate" : ""}</small></span>
        ${item.id === state.activeId ? `<span class="row-status" title="Active"><svg><use href="#i-check"></use></svg></span>` : ""}
      </button>`).join("")}
    </div>`).join("");
}

function updateShell() {
  document.querySelectorAll("[data-view]").forEach((button) => button.classList.toggle("is-active", button.dataset.view === state.view));
  document.querySelectorAll("[data-nav]").forEach((button) => {
    const isProjects = button.dataset.nav === "projects" && state.view === "projects";
    const isDesign = button.dataset.nav === "design" && state.view !== "projects";
    button.classList.toggle("is-active", isProjects || isDesign || button.dataset.nav === state.view);
    button.toggleAttribute("aria-current", isProjects || isDesign || button.dataset.nav === state.view);
  });
  document.querySelectorAll("[data-action=undo], [data-command=undo]").forEach((button) => { button.disabled = !state.undoStack.length; });
  document.querySelectorAll("[data-action=redo], [data-command=redo]").forEach((button) => { button.disabled = !state.redoStack.length; });
  const active = activeDocument();
  const documentName = document.querySelector("#active-document-name");
  if (documentName && active) documentName.textContent = active.name;
  if (!active) return;
  document.querySelectorAll("[data-document-name]").forEach((element) => { element.textContent = active.name; });
  document.querySelectorAll(".mobile-project-menu__title").forEach((element) => { element.textContent = active.name; });
  const subtitle = document.querySelector(".document-title p");
  if (subtitle) {
    subtitle.innerHTML = `<span>${new Intl.NumberFormat().format(active.sequence.length)} bp</span><i></i><span>${active.topology[0].toUpperCase()}${active.topology.slice(1)} DNA</span><i></i><span class="health-good"><b></b>Local revision r${active.revision}</span>`;
  }
  const documentTree = document.querySelector("[data-document-list]");
  if (documentTree) documentTree.innerHTML = renderDocumentTree();
  document.querySelectorAll("[data-view=features] em").forEach((badge) => { badge.textContent = active.features.length; });
  const documentCount = document.querySelector("#document-count");
  if (documentCount) documentCount.textContent = `${state.documents.length} design${state.documents.length === 1 ? "" : "s"}`;
}

function attachInteractiveMap() {
  const map = document.querySelector("#plasmid-map");
  if (!map) return;
  let start = null;
  let transform = { x: 0, y: 0, scale: 1 };
  const pointers = new Map();
  let pinch = null;
  const apply = () => { map.style.transform = `translate(${transform.x}px, ${transform.y}px) scale(${transform.scale})`; };
  const distance = () => {
    const [one, two] = [...pointers.values()];
    return one && two ? Math.hypot(one.x - two.x, one.y - two.y) : 0;
  };
  map.addEventListener("pointerdown", (event) => {
    if (event.target.closest("[data-feature-id], [data-site-position]")) return;
    pointers.set(event.pointerId, { x: event.clientX, y: event.clientY });
    if (pointers.size === 2) pinch = { distance: distance(), scale: transform.scale };
    else start = { x: event.clientX, y: event.clientY, transform: { ...transform } };
    map.setPointerCapture?.(event.pointerId);
  });
  map.addEventListener("pointermove", (event) => {
    if (!pointers.has(event.pointerId)) return;
    pointers.set(event.pointerId, { x: event.clientX, y: event.clientY });
    if (pointers.size >= 2 && pinch?.distance) {
      transform.scale = Math.max(.7, Math.min(1.8, pinch.scale * distance() / pinch.distance));
      apply();
      return;
    }
    if (!start) return;
    transform.x = start.transform.x + event.clientX - start.x;
    transform.y = start.transform.y + event.clientY - start.y;
    apply();
  });
  const releasePointer = (event) => {
    pointers.delete(event.pointerId);
    if (pointers.size === 1) {
      const remaining = [...pointers.values()][0];
      start = { x: remaining.x, y: remaining.y, transform: { ...transform } };
    } else if (!pointers.size) start = null;
    pinch = null;
  };
  map.addEventListener("pointerup", releasePointer);
  map.addEventListener("pointercancel", releasePointer);
  map.addEventListener("wheel", (event) => { event.preventDefault(); transform.scale = Math.max(.7, Math.min(1.8, transform.scale + (event.deltaY < 0 ? .08 : -.08))); apply(); }, { passive: false });
}

function attachSequenceTouchSelection() {
  const scroll = document.querySelector("#sequence-scroll");
  if (!scroll) return;
  let anchor = null;
  let dragged = false;
  const refreshBaseHighlights = () => {
    const active = activeDocument();
    scroll.querySelectorAll("[data-base-position]").forEach((base) => {
      base.classList.toggle("is-selected", isPositionSelected(Number(base.dataset.basePosition), active));
    });
  };
  scroll.addEventListener("pointerdown", (event) => {
    const base = event.target.closest("[data-base-position]");
    if (!base) return;
    anchor = Number(base.dataset.basePosition);
    dragged = false;
    scroll.setPointerCapture?.(event.pointerId);
  });
  scroll.addEventListener("pointermove", (event) => {
    if (anchor === null) return;
    const base = document.elementFromPoint(event.clientX, event.clientY)?.closest("[data-base-position]");
    if (!base) return;
    const position = Number(base.dataset.basePosition);
    if (position === anchor) return;
    dragged = true;
    state.selection = { start: Math.min(anchor, position), end: Math.max(anchor, position) + 1 };
    state.activeFeatureId = null;
    refreshBaseHighlights();
  });
  scroll.addEventListener("pointerup", () => {
    if (anchor === null) return;
    const didDrag = dragged;
    anchor = null;
    dragged = false;
    if (didDrag) render();
  });
  scroll.addEventListener("pointercancel", () => { anchor = null; dragged = false; });
}

function openModal(content) {
  const layer = document.querySelector("#modal-layer");
  if (!layer) return;
  layer.innerHTML = `<div class="modal-backdrop" data-action="close-modal"></div><section class="modal-sheet" role="dialog" aria-modal="true">${content}</section>`;
  document.body.classList.add("modal-open");
  requestAnimationFrame(() => layer.querySelector("input, textarea, select, button")?.focus());
}

function closeModal() {
  const layer = document.querySelector("#modal-layer");
  if (layer) layer.innerHTML = "";
  document.body.classList.remove("modal-open");
}

function openFeatureModal() {
  const document = activeDocument();
  const defaultStart = state.selection.start + 1;
  const defaultEnd = state.selection.end > state.selection.start ? state.selection.end : Math.min(document.sequence.length, defaultStart + 120);
  openModal(`<div class="modal-header"><div><span class="section-label">NEW ANNOTATION</span><h2>Add a feature</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><form id="feature-form" class="modal-form"><label>Name<input name="name" required placeholder="e.g. reporter CDS"/></label><div class="form-grid"><label>Type<select name="type"><option value="cds">CDS</option><option value="promoter">Promoter</option><option value="terminator">Terminator</option><option value="origin">Origin</option><option value="primer_bind">Primer binding site</option><option value="misc_feature">Misc. feature</option></select></label><label>Strand<select name="strand"><option value="1">Forward (+)</option><option value="-1">Reverse (−)</option></select></label></div><div class="form-grid"><label>Start <small>1-based</small><input type="number" name="start" min="1" max="${document.sequence.length}" value="${defaultStart}" required/></label><label>End <small>exclusive display</small><input type="number" name="end" min="1" max="${document.sequence.length}" value="${defaultEnd}" required/></label></div><label>Note<textarea name="note" rows="3" placeholder="Optional scientific note"></textarea></label><div class="modal-actions"><button type="button" class="soft-action" data-action="close-modal">Cancel</button><button type="submit" class="primary-action">Add annotation ${icon("arrow")}</button></div></form>`);
}

function openSequenceEditModal() {
  const document = activeDocument();
  const current = selectionSequence(document);
  const start = state.selection.start;
  const end = state.selection.end;
  openModal(`<div class="modal-header"><div><span class="section-label">SEQUENCE EDIT</span><h2>Edit bases</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><form id="sequence-edit-form" class="modal-form"><div class="edit-mode-picker"><button type="button" data-action="choose-edit-operation" data-operation="insert" class="${state.editOperation === "insert" ? "is-active" : ""}">Insert</button><button type="button" data-action="choose-edit-operation" data-operation="replace" class="${state.editOperation === "replace" ? "is-active" : ""}">Replace</button><button type="button" data-action="choose-edit-operation" data-operation="delete" class="${state.editOperation === "delete" ? "is-active" : ""}">Delete</button></div><div class="edit-target"><span>${icon("cursor")}</span><div><strong>${selectionLength(document) ? `${formatPosition(start)}–${formatPosition(end - 1)}` : `Cursor at ${formatPosition(start)}`}</strong><small>${selectionLength(document) ? `${formatBp(selectionLength(document))} currently selected` : "Tap or drag in Sequence view to select bases"}</small></div></div><label>DNA bases<textarea name="sequence" spellcheck="false" rows="5" placeholder="ACGT…" ${state.editOperation === "delete" ? "disabled" : ""}>${escapeHtml(state.editOperation === "insert" ? "" : current)}</textarea></label><p class="form-help">Only A, C, G, T, and N are retained. A feature’s coordinates are automatically remapped after the edit.</p><div class="modal-actions"><button type="button" class="soft-action" data-action="close-modal">Cancel</button><button type="submit" class="primary-action">Apply edit ${icon("check")}</button></div></form>`);
}

function openAssemblyModal(method = "gibson") {
  const document = activeDocument();
  const methodLabel = {
    gibson: "Gibson assembly", golden: "Golden Gate", ligation: "Restriction + ligation", synthetic: "Synthetic insert",
    gateway: "Gateway transfer", infusion: "In-Fusion / overlap assembly", ta: "TA / blunt cloning", oligo: "Oligo annealing",
  }[method] || "Assembly operation";
  const selection = selectionLength(document) ? `${formatPosition(state.selection.start)}–${formatPosition(state.selection.end - 1)}` : `cursor ${formatPosition(state.selection.start)}`;
  openModal(`<div class="modal-header"><div><span class="section-label">NEW PROVENANCE NODE</span><h2>${methodLabel}</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><form id="assembly-form" class="modal-form"><input type="hidden" name="method" value="${method}"/><div class="operation-input"><span>${icon("helix")}</span><div><strong>${escapeHtml(document.name)}</strong><small>${formatBp(document.sequence.length)} · ${document.topology} backbone · ${selection}</small></div><i>${icon("arrow")}</i></div><label>Operation name<input name="name" value="${methodLabel}" required/></label><label>Candidate product name<input name="productName" value="${escapeHtml(document.name)} · ${methodLabel}" required/></label><label>Optional insert / assembled sequence<textarea name="insertSequence" spellcheck="false" rows="4" placeholder="ACGT… (leave empty to record a plan only)"></textarea></label><label>Design note<textarea name="note" rows="3" placeholder="Describe fragments, enzymes, overlaps, or intended product"></textarea></label><p class="form-help">When an insert is supplied, this creates a separate draft product at the active selection/cursor. It is a sequence-planning aid, not a laboratory validation.</p><label class="checkbox-label"><input type="checkbox" name="confirm" required/> <span>I confirm the source construct remains unchanged and I will review warnings before laboratory use.</span></label><div class="modal-actions"><button type="button" class="soft-action" data-action="close-modal">Cancel</button><button type="submit" class="primary-action">Save plan / candidate ${icon("arrow")}</button></div></form>`);
}

function recordOperation(entry) {
  state.history.unshift({
    id: uid("event"),
    timestamp: new Date().toISOString(),
    warnings: [],
    ...entry,
  });
  state.history = state.history.slice(0, 24);
  persistLocalProject();
}

function getGoldenGateCuts(sequence, isCircular = false) {
  const searchable = isCircular ? sequence + sequence.slice(0, 10) : sequence;
  const cuts = [];
  const doubled = sequence + sequence;
  let index = searchable.indexOf("GGTCTC");
  while (index !== -1) {
    const t = (index + 7) % sequence.length;
    const overhang = doubled.slice(index + 7, index + 11);
    if (!cuts.some(c => c.t === t)) cuts.push({ t, overhang, type: "forward", pos: index % sequence.length });
    index = searchable.indexOf("GGTCTC", index + 1);
  }
  index = searchable.indexOf("GAGACC");
  while (index !== -1) {
    let t = (index - 5) % sequence.length;
    if (t < 0) t += sequence.length;
    const start = index - 5 < 0 ? index - 5 + sequence.length : index - 5;
    const overhang = doubled.slice(start, start + 4);
    if (!cuts.some(c => c.t === t)) cuts.push({ t, overhang, type: "reverse", pos: index % sequence.length });
    index = searchable.indexOf("GAGACC", index + 1);
  }
  return cuts.sort((a, b) => a.t - b.t);
}

function createAssemblyResult(form) {
  const source = activeDocument();
  const insert = normalizeSequence(form.insertSequence || "");
  const note = String(form.note || "").trim();
  if (!insert) {
    recordOperation({ type: "cloning", title: form.name, detail: note || `${form.method} planning record · no candidate product sequence supplied`, inputs: [source.id], output: null, warnings: ["Planning record only — no insert sequence was provided."] });
    closeModal();
    state.view = "history";
    render();
    showToast("Assembly plan saved in provenance");
    return;
  }
  const start = state.selection.start;
  const end = state.selection.end > state.selection.start ? state.selection.end : start;
  const product = clone(source);
  product.id = uid("doc");
  product.name = String(form.productName || `${source.name} candidate`).trim();
  product.description = `${form.name} candidate · derived locally from ${source.name}`;
  product.metadata = { ...product.metadata, parentDocumentId: source.id, operation: form.method, designNote: note, candidate: true };
  product.revision = 1;

  if (form.method === "golden") {
    const bbCuts = getGoldenGateCuts(source.sequence, source.topology === "circular");
    const testGG = (insSeq) => {
      const insCuts = getGoldenGateCuts(insSeq, false);
      if (bbCuts.length === 2 && insCuts.length === 2) {
        if (insCuts[0].overhang === bbCuts[0].overhang && insCuts[1].overhang === bbCuts[1].overhang) {
          return source.sequence.slice(0, bbCuts[0].t) + insSeq.slice(insCuts[0].t, insCuts[1].t) + source.sequence.slice(bbCuts[1].t);
        }
        if (insCuts[0].overhang === bbCuts[1].overhang && insCuts[1].overhang === bbCuts[0].overhang) {
          return source.sequence.slice(bbCuts[0].t, bbCuts[1].t) + insSeq.slice(insCuts[0].t, insCuts[1].t);
        }
      }
      return null;
    };

    let assembledSeq = testGG(insert);
    if (!assembledSeq) assembledSeq = testGG(reverseComplement(insert));

    if (!assembledSeq) {
      showToast("Golden Gate failed: incompatible BsaI overhangs or incorrect site counts.");
      return;
    }
    product.sequence = assembledSeq;
    product.topology = "circular";
    product.features = [];
    state.selection = { start: 0, end: 0 };
  } else {
    product.sequence = `${source.sequence.slice(0, start)}${insert}${source.sequence.slice(end)}`;
    remapAfterEdit(product, start, end - start, insert.length);
    state.selection = { start, end: start + insert.length };
  }

  state.documents.push(product);
  state.activeId = product.id;
  state.activeFeatureId = null;
  const warnings = [
    "Candidate sequence created from the supplied insert; validate junctions, enzyme compatibility, and laboratory parameters.",
  ];
  recordOperation({ type: "cloning", title: form.name, detail: `${source.name} + ${formatBp(insert.length)} insert → ${product.name}`, inputs: [source.id], output: product.id, parameters: { method: form.method, insertionStart: start, replacedLength: end - start, insertLength: insert.length, sourceHash: simpleHash(source.sequence), productHash: simpleHash(product.sequence) }, warnings });
  closeModal();
  state.view = "history";
  render();
  showToast(`Candidate ${product.name} created`);
}

function openReadModal() {
  openModal(`<div class="modal-header"><div><span class="section-label">ADD LOCAL READ</span><h2>Paste a sequence read</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><form id="read-form" class="modal-form"><label>Read name<input name="name" required placeholder="e.g. sample-01 forward read"/></label><label>Sequence<textarea name="sequence" spellcheck="false" rows="6" required placeholder="ACGT…"></textarea></label><label>Source / note<input name="source" placeholder="e.g. manually pasted, FASTA read"/></label><p class="form-help">This portable workspace accepts nucleotide reads. The Python sequencing package offers an optional real-AB1 parser boundary; this browser view never fabricates chromatograms.</p><div class="modal-actions"><button type="button" class="soft-action" data-action="close-modal">Cancel</button><button type="submit" class="primary-action">Add read ${icon("arrow")}</button></div></form>`);
}

function openCodonOptimizeModal() {
  openModal(`<div class="modal-header"><div><span class="section-label">CODING SEQUENCE DESIGN</span><h2>Reverse translate protein</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><form id="codon-optimize-form" class="modal-form"><label>Candidate name<input name="name" required value="Synthetic coding candidate"/></label><label>Protein sequence<textarea name="protein" spellcheck="false" rows="5" required placeholder="M... (one-letter amino acids; * is a stop codon)"></textarea></label><label>Design note<textarea name="note" rows="3" placeholder="Optional intended host, assay, or ordering note"></textarea></label><p class="form-help">Uses a deterministic generic codon preference table. It does not provide organism-specific optimisation, expression guarantees, safety screening, or ordering approval.</p><label class="checkbox-label"><input type="checkbox" name="confirm" required/> <span>I will independently review the generated sequence and intended use before any ordering or laboratory work.</span></label><div class="modal-actions"><button type="button" class="soft-action" data-action="close-modal">Cancel</button><button type="submit" class="primary-action">Create DNA candidate ${icon("arrow")}</button></div></form>`);
}

function createCodonCandidate(values) {
  const protein = String(values.protein || "").toUpperCase().replace(/\s+/g, "");
  if (!protein || /[^ACDEFGHIKLMNPQRSTVWY*]/.test(protein)) return showToast("Use standard one-letter amino acids (and optional * stop codon)");
  const sequence = [...protein].map((aminoAcid) => PREFERRED_CODON[aminoAcid]).join("");
  if (!sequence) return showToast("No DNA codons could be generated");
  const source = activeDocument();
  const name = String(values.name || "Synthetic coding candidate").trim() || "Synthetic coding candidate";
  const candidate = {
    id: uid("doc"),
    name,
    description: `Generic reverse-translation candidate · ${protein.length} amino acids`,
    sequence,
    topology: "linear",
    moleculeType: "DNA",
    revision: 1,
    metadata: {
      candidate: true,
      parentDocumentId: source?.id || null,
      operation: "generic_reverse_translation",
      designNote: String(values.note || "").trim(),
      proteinSequence: protein,
      codonTable: "browser-generic-preference-v1",
    },
    features: [feature("cds", name, 0, sequence.length, 1, {
      translation: protein.replace(/\*$/, ""),
      note: "Generic reverse-translation candidate; organism-specific optimisation and review required.",
    })],
  };
  state.documents.push(candidate);
  state.activeId = candidate.id;
  state.selection = { start: 0, end: sequence.length };
  state.activeFeatureId = candidate.features[0].id;
  recordOperation({
    type: "design",
    title: "Created reverse-translation candidate",
    detail: `${protein.length} amino acids → ${formatBp(sequence.length)} DNA using a generic preference table`,
    inputs: source ? [source.id] : [],
    output: candidate.id,
    parameters: { proteinLength: protein.length, codonTable: "browser-generic-preference-v1" },
    warnings: ["Generic codon usage only — select organism-specific preferences and independently review before use."],
  });
  closeModal();
  state.view = "history";
  render();
  showToast(`Created ${candidate.name}`);
}

function addRead(values) {
  const sequence = normalizeSequence(values.sequence);
  if (sequence.length < 12) return showToast("A sequence read needs at least 12 valid DNA bases");
  state.reads.push({ id: uid("read"), name: String(values.name).trim() || "Untitled read", sequence, source: String(values.source || "manual read").trim(), quality: Array.from({ length: sequence.length }, () => 30), createdAt: new Date().toISOString() });
  recordOperation({ type: "analysis", title: "Added verification read", detail: `${values.name || "Untitled read"} · ${formatBp(sequence.length)}`, inputs: [activeDocument().id], output: null });
  closeModal();
  state.activeTool = "verification";
  state.view = "tools";
  render();
  showToast("Read added locally");
}

function createConsensusDocument() {
  const source = activeDocument();
  const summary = verificationSummary(source);
  if (!state.reads.length) return showToast("Add at least one read first");
  const product = clone(source);
  product.id = uid("doc");
  product.name = `${source.name} consensus`;
  product.description = `Consensus candidate from ${state.reads.length} locally aligned read${state.reads.length === 1 ? "" : "s"}`;
  product.sequence = summary.consensus;
  product.revision = 1;
  product.metadata = { ...source.metadata, parentDocumentId: source.id, consensus: true, reads: state.reads.map((read) => read.id) };
  state.documents.push(product);
  state.activeId = product.id;
  state.verificationRun = { sourceId: source.id, outputId: product.id, ranAt: new Date().toISOString(), variants: summary.variants, meanIdentity: summary.meanIdentity };
  recordOperation({ type: "verification", title: "Created consensus candidate", detail: `${state.reads.length} reads · ${(summary.meanIdentity * 100).toFixed(1)}% mean identity · ${summary.variants.length} confident substitutions`, inputs: [source.id, ...state.reads.map((read) => read.id)], output: product.id, warnings: ["Consensus is a local algorithmic candidate. Review trace quality and alignments before accepting a sequence."] });
  state.view = "history";
  render();
  showToast(`Consensus candidate ${product.name} created`);
}

function runVerification() {
  const document = activeDocument();
  if (!state.reads.length) return openReadModal();
  const summary = verificationSummary(document);
  state.verificationRun = { sourceId: document.id, ranAt: new Date().toISOString(), variants: summary.variants, meanIdentity: summary.meanIdentity, coveredBases: summary.coveredBases };
  recordOperation({ type: "verification", title: "Aligned sequence reads", detail: `${state.reads.length} reads · ${(summary.meanIdentity * 100).toFixed(1)}% mean identity · ${(summary.coveredBases / document.sequence.length * 100).toFixed(1)}% reference coverage`, inputs: [document.id, ...state.reads.map((read) => read.id)], output: null, warnings: summary.variants.length ? ["Confident substitutions need read/trace review before sequence acceptance."] : [] });
  render();
  showToast("Read alignment recorded");
}

function focusVariant(position) {
  const document = activeDocument();
  state.selection = { start: position, end: Math.min(document.sequence.length, position + 1) };
  state.activeFeatureId = null;
  state.view = "sequence";
  render();
  setTimeout(scrollToSelection, 80);
}

function inspectRead(id) {
  const read = state.reads.find((candidate) => candidate.id === id);
  if (!read) return;
  const alignment = bestUngappedAlignment(activeDocument().sequence, read.sequence, activeDocument().topology === "circular");
  openModal(`<div class="modal-header"><div><span class="section-label">READ INSPECTION</span><h2>${escapeHtml(read.name)}</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><div class="history-detail"><dl class="feature-dl"><div><dt>Length</dt><dd>${formatBp(read.sequence.length)}</dd></div><div><dt>Mapping</dt><dd>${alignment ? `${formatPosition(alignment.start)} · ${alignment.strand} strand` : "unmapped"}</dd></div><div><dt>Identity</dt><dd>${alignment ? `${(alignment.identity * 100).toFixed(2)}%` : "—"}</dd></div><div><dt>Substitutions</dt><dd>${alignment?.variants.length || 0}</dd></div></dl><code class="read-sequence">${escapeHtml(read.sequence)}</code>${alignment?.variants.length ? `<div class="warning-box">${icon("info")} ${alignment.variants.slice(0, 5).map((variant) => `${formatPosition(variant.position)} ${variant.reference}→${variant.alternate}`).join(" · ")}</div>` : ""}<div class="modal-actions"><button class="primary-action" data-action="close-modal">Done</button></div></div>`);
}

function showToast(message) {
  const toast = document.querySelector("#toast");
  if (!toast) return;
  toast.textContent = message;
  toast.classList.add("is-visible");
  clearTimeout(showToast.timer);
  showToast.timer = setTimeout(() => toast.classList.remove("is-visible"), 2800);
}

function handleFile(file) {
  if (!file) return;
  if (file.name.match(/\.ogs$/i)) {
    showToast("This browser shell opens .ogs.json snapshots. Use the OpenGene CLI for ZIP + SQLite .ogs projects.");
    return;
  }
  const reader = new FileReader();
  reader.onload = () => {
    try {
      const text = String(reader.result || "");
      let document;
      if (file.name.match(/\.(fastq|fq)$/i) || /^@[^\n]+\n[^\n]+\n\+[^\n]*\n/m.test(text)) {
        const reads = parseFastq(text, file.name);
        state.reads.push(...reads);
        recordOperation({ type: "import", title: `Imported ${reads.length} sequence read${reads.length === 1 ? "" : "s"}`, detail: `FASTQ reads from ${file.name}`, inputs: [activeDocument().id], output: null });
        state.activeTool = "verification";
        state.view = "tools";
        render();
        showToast(`Imported ${reads.length} read${reads.length === 1 ? "" : "s"}`);
        return;
      }
      if (file.name.match(/\.(gb|gbk|genbank)$/i) || /^LOCUS\s+/m.test(text)) document = parseGenBank(text, file.name);
      else if (file.name.match(/\.json$/i)) {
        const project = JSON.parse(text);
        if (!Array.isArray(project.documents)) throw new Error("This is not an OpenGene project snapshot.");
        state.documents = project.documents;
        state.activeId = project.activeId || project.documents[0]?.id;
        state.history = project.history || [];
        state.reads = Array.isArray(project.reads) ? project.reads : [];
        state.verificationRun = project.verificationRun || null;
        state.undoStack = [];
        state.redoStack = [];
        state.selection = { start: 0, end: 0 };
        persistLocalProject();
        render();
        showToast(`Opened ${project.documents.length} local design${project.documents.length === 1 ? "" : "s"}`);
        return;
      } else document = parseFasta(text, file.name);
      state.documents.push(document);
      state.activeId = document.id;
      state.selection = { start: 0, end: 0 };
      state.activeFeatureId = null;
      recordHistory({ type: "import", title: `Imported ${document.name}`, detail: `${formatBp(document.sequence.length)} from ${file.name}` });
      render();
      showToast(`Imported ${document.name}`);
    } catch (error) {
      showToast(error.message || "Could not import that file");
    }
  };
  reader.readAsText(file);
}

function eventTargetAction(target) {
  return target.closest("[data-action], [data-command], [data-tool], [data-view], [data-view-target], [data-feature-id], [data-base-position], [data-site-position]");
}

function attachEvents() {
  document.addEventListener("click", (event) => {
    const control = eventTargetAction(event.target);
    if (!control) return;
    const activeDoc = activeDocument();
    if (control.dataset.viewTarget) {
      const target = control.dataset.viewTarget.replace("-view", "");
      state.view = target === "files" ? "projects" : target === "tools" ? "tools" : target;
      render();
      return;
    }
    if (control.dataset.tool && !control.dataset.action) {
      state.activeTool = ({ "restriction-digest": "restriction", pcr: "primers", gibson: "cloning", align: "verification" })[control.dataset.tool] || "restriction";
      state.view = "tools";
      render();
      return;
    }
    if (control.dataset.command) {
      const commandAction = {
        import: "import", export: "export-genbank", save: "save-project", undo: "undo", redo: "redo",
        annotate: "add-feature", clone: "create-assembly", "find-sites": "open-restriction",
        "new-document": "import", "open-sidebar": "toggle-sidebar", "focus-feature": "zoom-selection",
      }[control.dataset.command];
      if (commandAction) control.dataset.action = commandAction;
      else if (control.dataset.command === "search") { state.view = "features"; render(); setTimeout(() => window.document.querySelector("#feature-search")?.focus(), 20); return; }
      else if (control.dataset.command === "share") { saveProject(); return; }
      else return;
    }
    if (control.dataset.view) {
      state.view = control.dataset.view;
      render();
      return;
    }
    if (control.dataset.featureId && !control.dataset.action) { selectFeature(control.dataset.featureId); return; }
    if (control.dataset.basePosition !== undefined) {
      const position = Number(control.dataset.basePosition);
      if (event.shiftKey) state.selection.end = position + 1;
      else state.selection = { start: position, end: position };
      const containing = activeDoc.features.find((item) => featureContains(item, position, activeDoc));
      state.activeFeatureId = containing?.id || null;
      render();
      return;
    }
    if (control.dataset.sitePosition !== undefined) {
      const position = Number(control.dataset.sitePosition);
      state.selection = { start: position, end: Math.min(activeDoc.sequence.length, position + 6) };
      showToast(`${control.dataset.siteName} at ${formatPosition(position)}`);
      render();
      return;
    }
    const action = control.dataset.action;
    switch (action) {
      case "import": window.document.getElementById("file-input")?.click(); break;
      case "toggle-sidebar": document.body.classList.toggle("sidebar-open"); break;
      case "add-feature": openFeatureModal(); break;
      case "edit-sequence": openSequenceEditModal(); break;
      case "close-modal": closeModal(); break;
      case "undo": undo(); break;
      case "redo": redo(); break;
      case "toggle-sites": state.showSites = !state.showSites; render(); break;
      case "set-map-mode": state.mode = control.dataset.mode; render(); break;
      case "toggle-map-theme": state.darkMap = !state.darkMap; render(); break;
      case "map-zoom-in": zoomMap(.12); break;
      case "map-zoom-out": zoomMap(-.12); break;
      case "map-reset": resetMap(); break;
      case "open-tool": state.activeTool = control.dataset.tool; state.view = "tools"; render(); break;
      case "open-document": openDocument(control.dataset.documentId); break;
      case "duplicate-document": duplicateActiveDocument(); break;
      case "delete-feature": deleteFeature(control.dataset.featureId); break;
      case "feature-menu": selectFeature(control.dataset.featureId); showToast("Feature selected — use the inspector to copy or focus it"); break;
      case "focus-feature": selectFeature(control.dataset.featureId); state.view = "map"; render(); break;
      case "copy-feature-sequence": copyFeatureSequence(control.dataset.featureId); break;
      case "copy-selection": copyText(selectionSequence(activeDoc), selectionLength(activeDoc) ? "Selected sequence copied" : "Nothing selected"); break;
      case "zoom-selection": state.view = "sequence"; render(); setTimeout(scrollToSelection, 100); break;
      case "reverse-complement-selection": reverseComplementSelection(); break;
      case "translate-selection": translateSelection(); break;
      case "design-primers": showToast("Primer suggestions refreshed around the current selection"); render(); break;
      case "copy-primer": copyText(control.dataset.primer, "Primer copied"); break;
      case "annotate-primers": annotatePrimers(); break;
      case "add-orf": addBrowserOrf(Number(control.dataset.orfIndex)); break;
      case "annotate-all-orfs": annotateBrowserOrfs(); break;
      case "auto-annotate": autoAnnotate(); break;
      case "open-codon-optimize": openCodonOptimizeModal(); break;
      case "add-read": openReadModal(); break;
      case "run-verification": runVerification(); break;
      case "create-consensus": createConsensusDocument(); break;
      case "inspect-read": inspectRead(control.dataset.readId); break;
      case "focus-variant": focusVariant(Number(control.dataset.position)); break;
      case "copy-translation": copyText(control.dataset.protein, "Protein translation copied"); break;
      case "enzyme-toggle": showToast(`${control.dataset.enzyme} is pinned for this preview`); break;
      case "run-digest": recordHistory({ type: "analysis", title: "Ran restriction digest", detail: "EcoRI + BamHI + HindIII local analysis" }); showToast("Digest analysis recorded"); render(); break;
      case "select-fragment": state.selection = { start: Number(control.dataset.start), end: Number(control.dataset.end) }; state.view = "sequence"; render(); setTimeout(scrollToSelection, 100); break;
      case "create-assembly": openAssemblyModal(control.dataset.method || "gibson"); break;
      case "run-gel": recordHistory({ type: "analysis", title: "Rendered virtual gel", detail: "Deterministic fragment-size visualization" }); showToast("Virtual gel rendered locally"); render(); break;
      case "export-fasta": exportFasta(); break;
      case "export-genbank": exportGenBank(); break;
      case "save-project": saveProject(); break;
      case "export-manifest": download(`${safeFileName(activeDoc.name)}-manifest.json`, JSON.stringify(projectSnapshot(), null, 2), "application/json"); showToast("Project manifest downloaded"); break;
      case "open-restriction": state.activeTool = "restriction"; state.view = "tools"; render(); break;
      case "history-detail": historyDetail(control.dataset.historyId); break;
      case "choose-edit-operation": state.editOperation = control.dataset.operation; openSequenceEditModal(); break;
      default: break;
    }
  });

  document.addEventListener("input", (event) => {
    if (event.target.id === "feature-search") {
      state.search = event.target.value;
      const body = document.querySelector("#feature-table-body");
      if (body) {
        const active = activeDocument();
        body.innerHTML = renderFeatureTableRows(active, active.features.find((item) => item.id === state.activeFeatureId));
      }
    }
  });

  document.addEventListener("submit", (event) => {
    if (event.target.id === "feature-form") {
      event.preventDefault();
      addFeatureFromForm(Object.fromEntries(new FormData(event.target).entries()));
    }
    if (event.target.id === "sequence-edit-form") {
      event.preventDefault();
      const form = Object.fromEntries(new FormData(event.target).entries());
      const document = activeDocument();
      const start = state.selection.start;
      const end = state.selection.end;
      const operation = state.editOperation;
      if (operation === "insert") replaceSequenceRegion(start, start, form.sequence || "");
      else if (operation === "delete") replaceSequenceRegion(start, end > start ? end : Math.min(document.sequence.length, start + 1), "");
      else replaceSequenceRegion(start, end > start ? end : Math.min(document.sequence.length, start + 1), form.sequence || "");
      closeModal();
    }
    if (event.target.id === "assembly-form") {
      event.preventDefault();
      const form = Object.fromEntries(new FormData(event.target).entries());
      createAssemblyResult(form);
    }
    if (event.target.id === "read-form") {
      event.preventDefault();
      addRead(Object.fromEntries(new FormData(event.target).entries()));
    }
    if (event.target.id === "codon-optimize-form") {
      event.preventDefault();
      createCodonCandidate(Object.fromEntries(new FormData(event.target).entries()));
    }
  });

  document.addEventListener("keydown", (event) => {
    if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === "z") { event.preventDefault(); event.shiftKey ? redo() : undo(); }
    if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === "k") { event.preventDefault(); state.view = "features"; render(); setTimeout(() => document.querySelector("#feature-search")?.focus(), 20); }
    if (event.key === "Escape") closeModal();
  });

  document.getElementById("file-input")?.addEventListener("change", (event) => {
    handleFile(event.target.files?.[0]);
    event.target.value = "";
  });
}

function zoomMap(delta) {
  const map = document.querySelector("#plasmid-map");
  if (!map) return;
  const current = Number(map.dataset.scale || 1);
  const next = Math.max(.75, Math.min(1.8, current + delta));
  map.dataset.scale = next;
  map.style.transform = `scale(${next})`;
}

function resetMap() {
  const map = document.querySelector("#plasmid-map");
  if (map) { map.dataset.scale = 1; map.style.transform = ""; }
}

function scrollToSelection() {
  const scroll = document.querySelector("#sequence-scroll");
  const row = document.querySelector(`[data-sequence-row="${Math.floor(state.selection.start / 60) * 60}"]`);
  if (scroll && row) row.scrollIntoView({ behavior: "smooth", block: "center" });
}

function copyText(value, message) {
  if (!value) return showToast(message || "Nothing to copy");
  navigator.clipboard?.writeText(value).then(() => showToast(message || "Copied to clipboard")).catch(() => showToast(value.length > 72 ? "Selection ready to copy" : value));
}

function copyFeatureSequence(id) {
  const document = activeDocument();
  const item = document.features.find((candidate) => candidate.id === id);
  if (!item) return;
  const value = featureSegments(item, document.sequence.length).map(([start, end]) => document.sequence.slice(start, end)).join("");
  copyText(item.strand === -1 ? reverseComplement(value) : value, `${item.name} sequence copied`);
}

function reverseComplementSelection() {
  const document = activeDocument();
  const value = selectionSequence(document);
  if (!value || state.selection.end < state.selection.start) return showToast("Select a non-wrapping range before reverse complementing it");
  replaceSequenceRegion(state.selection.start, state.selection.end, reverseComplement(value));
}

function translateSelection() {
  const sequence = selectionSequence(activeDocument());
  if (sequence.length < 3) return showToast("Select at least one codon to translate");
  const code = GENETIC_CODE;
  let protein = "";
  for (let index = 0; index + 2 < sequence.length; index += 3) protein += code[sequence.slice(index, index + 3)] || "X";
  openModal(`<div class="modal-header"><div><span class="section-label">TRANSLATION</span><h2>${Math.floor(sequence.length / 3)} codons</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><div class="translation-result"><code>${escapeHtml(protein)}</code><p>Translated in frame 1 from the active forward-strand selection. Stop codons are shown as <strong>*</strong>.</p><div class="modal-actions"><button class="soft-action" data-action="close-modal">Close</button><button class="primary-action" data-action="copy-translation" data-protein="${protein}">${icon("copy")} Copy protein</button></div></div>`);
}

function annotatePrimers() {
  const document = activeDocument();
  const candidates = primerCandidates(document);
  commit("Added primer annotations", () => {
    const [forward, reverse] = candidates;
    document.features.push(feature("primer_bind", forward.name, forward.position, forward.position + forward.sequence.length, 1, { note: `Tm ${forward.tm}°C · GC ${forward.gc}%` }));
    document.features.push(feature("primer_bind", reverse.name, reverse.position, reverse.position + reverse.sequence.length, -1, { note: `Tm ${reverse.tm}°C · GC ${reverse.gc}%` }));
  }, "Primer pair from local heuristic");
}

function exportFasta() {
  const document = activeDocument();
  download(`${safeFileName(document.name)}.fasta`, fastaText(document), "text/x-fasta");
  showToast("FASTA downloaded");
}

function exportGenBank() {
  const document = activeDocument();
  download(`${safeFileName(document.name)}.gb`, genbankText(document), "text/plain");
  showToast("GenBank downloaded");
}

function saveProject() {
  const document = activeDocument();
  persistLocalProject();
  download(`${safeFileName(document?.name || "opengene-project")}.ogs.json`, JSON.stringify(projectSnapshot(), null, 2), "application/json");
  showToast("Local project snapshot downloaded");
}

function safeFileName(value) {
  return String(value || "opengene").toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/(^-|-$)/g, "") || "opengene";
}

function historyDetail(id) {
  const entry = state.history.find((candidate) => candidate.id === id);
  if (!entry) return;
  openModal(`<div class="modal-header"><div><span class="section-label">PROVENANCE NODE</span><h2>${escapeHtml(entry.title)}</h2></div><button class="icon-action" data-action="close-modal" aria-label="Close">${icon("close")}</button></div><div class="history-detail"><dl class="feature-dl"><div><dt>Recorded</dt><dd>${new Date(entry.timestamp).toLocaleString()}</dd></div><div><dt>Operation type</dt><dd>${escapeHtml(entry.type)}</dd></div><div><dt>Input documents</dt><dd>${entry.inputs.length || "None"}</dd></div><div><dt>Output document</dt><dd>${entry.output ? "Active project construct" : "No product"}</dd></div><div><dt>Parameters</dt><dd>Captured in local record</dd></div></dl><p>${escapeHtml(entry.detail || "No additional detail.")}</p>${entry.warnings?.length ? `<div class="warning-box">${icon("info")} ${escapeHtml(entry.warnings.join(" "))}</div>` : ""}<div class="modal-actions"><button class="primary-action" data-action="close-modal">Done</button></div></div>`);
}

function icon(name) {
  const paths = {
    plus: '<path d="M12 5v14M5 12h14"/>',
    edit: '<path d="m4 16.5-.8 3.3 3.3-.8L17.4 8.1l-2.6-2.6L4 16.5Z"/><path d="m13.5 6.8 2.6 2.6"/>',
    scissors: '<path d="M6.5 8.5 17.8 15M6.5 15l11.3-6.5"/><circle cx="5.5" cy="7" r="2.2"/><circle cx="5.5" cy="17" r="2.2"/>',
    circle: '<circle cx="12" cy="12" r="7"/>',
    line: '<path d="M4 12h16"/><circle cx="5" cy="12" r="1"/><circle cx="19" cy="12" r="1"/>',
    focus: '<path d="M8 4H4v4M16 4h4v4M20 16v4h-4M4 16v4h4"/><circle cx="12" cy="12" r="4"/>',
    sun: '<circle cx="12" cy="12" r="3.5"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/>',
    hand: '<path d="M8 12V6.8a1.3 1.3 0 0 1 2.6 0V11M10.6 11V5.6a1.3 1.3 0 0 1 2.6 0V11M13.2 11V6.8a1.3 1.3 0 0 1 2.6 0V13M15.8 12v-2a1.3 1.3 0 0 1 2.6 0v5.2c0 3.5-2.8 5.8-6.2 5.8H11c-1.8 0-3.3-1.1-4.2-2.5L4.6 15.3a1.4 1.4 0 0 1 2.1-1.8L8 15"/>',
    arrow: '<path d="M5 12h14M14 6l6 6-6 6"/>',
    tag: '<path d="M3.5 12.3V5h7.3l8.7 8.7-6.5 6.5-8.8-8.8Z"/><circle cx="8" cy="8.3" r="1"/>',
    cursor: '<path d="m5 3 12 8-5 1.5L10.5 18 5 3Z"/><path d="m13 14 4 5"/>',
    copy: '<rect x="8" y="8" width="11" height="12" rx="1.5"/><path d="M16 8V5.5A1.5 1.5 0 0 0 14.5 4h-9A1.5 1.5 0 0 0 4 5.5v10A1.5 1.5 0 0 0 5.5 17H8"/>',
    reverse: '<path d="M7 7H4l3-3M4 7c1.7-2 4.3-3 7-3 4.4 0 8 3.6 8 8M17 17h3l-3 3M20 17c-1.7 2-4.3 3-7 3-4.4 0-8-3.6-8-8"/>',
    translate: '<path d="M4 5h10M9 5c0 6-2 10-5 12M6 11c1.2 2 2.7 3.8 4.5 5"/><path d="M15 19l3-10 3 10M16 16h4"/>',
    trash: '<path d="M5 7h14M10 11v5M14 11v5M8 7l.8-3h6.4l.8 3M7 7l.8 13h8.4L17 7"/>',
    search: '<circle cx="10.5" cy="10.5" r="5.5"/><path d="m15 15 4.5 4.5"/>',
    dots: '<circle cx="5" cy="12" r="1"/><circle cx="12" cy="12" r="1"/><circle cx="19" cy="12" r="1"/>',
    primer: '<path d="M5 7h14M5 17h14M7 5v4M11 5v4M15 15v4M19 15v4"/>',
    assembly: '<circle cx="6" cy="12" r="2.5"/><circle cx="18" cy="6" r="2.5"/><circle cx="18" cy="18" r="2.5"/><path d="m8.2 10.8 7.5-3.6M8.2 13.2l7.5 3.6"/>',
    gel: '<rect x="5" y="3" width="14" height="18" rx="2"/><path d="M8 8h8M8 12h8M8 16h8"/>',
    verify: '<path d="M5 4h11l3 3v13H5V4Z"/><path d="M8 11h5M8 15h3M15 14l1.5 1.5L20 12"/>',
    spark: '<path d="m12 3 1.5 5.5L19 10l-5.5 1.5L12 17l-1.5-5.5L5 10l5.5-1.5L12 3ZM19 15l.7 2.3L22 18l-2.3.7L19 21l-.7-2.3L16 18l2.3-.7L19 15Z"/>',
    play: '<path d="m8 5 10 7-10 7V5Z"/>',
    check: '<path d="m5 12 4.2 4.2L19 6.5"/>',
    history: '<path d="M4 12a8 8 0 1 0 2.3-5.6L4 8.7M4 4v4.7h4.7M12 7v5l3.5 2"/>',
    upload: '<path d="M12 16V4M7.5 8.5 12 4l4.5 4.5M5 16v3h14v-3"/>',
    download: '<path d="M12 4v12M7.5 11.5 12 16l4.5-4.5M5 19v1h14v-1"/>',
    pin: '<path d="m12 21 5-7a6 6 0 1 0-10 0l5 7Z"/><circle cx="12" cy="10" r="2"/>',
    info: '<circle cx="12" cy="12" r="8"/><path d="M12 11v5M12 8h.01"/>',
    helix: '<path d="M6 4c8 0 4 16 12 16M18 4C10 4 14 20 6 20M7 7h10M7 17h10"/>',
    close: '<path d="m6 6 12 12M18 6 6 18"/>',
  };
  return `<svg class="icon" viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">${paths[name] || paths.circle}</svg>`;
}

function bootstrap() {
  if (!restoreLocalProject()) {
    const demo = createDemoDocument();
    state.documents = [demo];
    state.activeId = demo.id;
    state.reads = createDemoReads(demo);
    createInitialHistory();
  }
  attachEvents();
  render();
  if ("serviceWorker" in navigator) navigator.serviceWorker.register("./service-worker.js").catch(() => {});
}

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", bootstrap);
else bootstrap();
