USE kinesin;

-- changed preferences out of 'safe updates' mode to allow deleting all contents; may want to put back (prefs, sql editor, reconnect to db) when adding real data.

 -- DELETE FROM kinesin_family;
 -- DELETE FROM gene;
DELETE FROM mutation;
DELETE FROM source_info;
DELETE FROM impact;
DELETE FROM tissue;
