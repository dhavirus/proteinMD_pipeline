// Unsent work (decisions, regions, variants) kept in this browser only, per input file.
// It is a convenience, not the record: the exported manifest is the record.

const key = (inputSha256) => `simprep-draft-v1-${inputSha256}`;

export function loadDraft(inputSha256) {
  try {
    const raw = localStorage.getItem(key(inputSha256));
    return raw ? JSON.parse(raw) : null;
  } catch {
    return null;
  }
}

export function saveDraft(inputSha256, state, savedAt) {
  try {
    localStorage.setItem(key(inputSha256), JSON.stringify({ savedAt, decisions: state.decisions, regions: state.regions, variants: state.variants }));
  } catch {
    // Storage unavailable (private window, blocked site data): the page still works.
  }
}

export function clearDraft(inputSha256) {
  try {
    localStorage.removeItem(key(inputSha256));
  } catch {
    // Nothing to clear.
  }
}
